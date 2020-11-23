#include <iostream>
#include <vector>
#include <cmath>

#include "gdal_priv.h"
#include "cpl_conv.h"

#include "lasreader.hpp"


struct Coordinate {
    double id, x, y;

    Coordinate(int paramID, double paramX, double paramY) : id(paramID), x(paramX), y(paramY) {}
};

struct Timings {
    int firstTime, lastTime;

    Timings(int first, int last) : firstTime(first), lastTime(last) {}
};

struct Bbox {
    int minX, minY, maxX, maxY;

    int xDiff = maxX - minX;
    int yDiff = maxY - minY;

    Bbox(int min_x, int min_y, int max_x, int max_y) : minX(min_x), minY(min_y), maxX(max_x), maxY(max_y) {}
};

int roundUp(F64 number) {
    int returnValue = (int) number;
    if (returnValue % 10 == 9){ // Last digit is a 9
        returnValue++;
    }
    return returnValue;
}

int main(int argc, char **argv) {
    if (argc != 5) {
        std::cout << "Invalid number of input arguments!\n";
        std::cout << "arg1: input_file, arg2: output_file, arg3: cell_count, arg4: thinning_factor\n";
        return 0;
    }

    const char *inputFile = argv[1];
    const char *outputFile = argv[2];
    const int CELL_COUNT = std::stoi(argv[3]); // std::stoi = cast to int
    const int THINNING_FACTOR = std::stoi(argv[4]);

    std::vector<Coordinate> leftTopCorners;
    std::vector<Timings> timings;

    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(inputFile);
    LASreader *lasreader = lasreadopener.open();

    // In AHN3 all corner points are integers, also round them up if precision is off (84999 -> 85000)
    const Bbox bbox = Bbox(roundUp(lasreader->get_min_x()), roundUp(lasreader->get_min_y()),
                           roundUp(lasreader->get_max_x()), roundUp(lasreader->get_max_y()));

    std::cout << "BBox boundaries: minX = " << bbox.minX << ", minY = " << bbox.minY << ", maxX = " << bbox.maxX
              << ", maxY = " << bbox.maxY << "\n";

    const int numPoints = lasreader->npoints;

    std::cout << "Number of points: " << numPoints << "\n";

    const int xCellWidth = bbox.xDiff / CELL_COUNT;
    const int yCellWidth = bbox.yDiff / CELL_COUNT;

    int count = 0;

    for (int x = bbox.minX; x < bbox.maxX; x += xCellWidth) {
        for (int y = bbox.minY; y < bbox.maxY; y += yCellWidth) {
            leftTopCorners.emplace_back(count, x, y);
            count++;
        }
    }

    int streamTime = 0;
    int lastPercentage = -1;

    while (lasreader->read_point()) {

        streamTime++;

        if (streamTime % 20000 == 0) {
            double percentage = round((double) streamTime / (double) numPoints * 100);
            if (percentage != lastPercentage) {
                std::cout << percentage << "% done!\n";
                lastPercentage = percentage;
            }
        }

        if (streamTime % THINNING_FACTOR == 0) {

            for (Coordinate corner : leftTopCorners) {

                if (lasreader->point.inside_rectangle(corner.x, corner.y, corner.x + xCellWidth, corner.y + yCellWidth)) {

                    try { // Vector exists, update last entry time
                        timings.at(corner.id).lastTime = streamTime;
                    }
                    catch (...) { // Vector doesn't exist, set first + last time and push back
                        timings.emplace_back(streamTime, streamTime);
                    }

                    break;
                }
            }
        }

    }

    std::cout << "Writing GeoTIFF \n";

    GDALAllRegister();

    const char *pszFormat = "GTiff";
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

    if (poDriver == nullptr)
        exit(1);

    GDALDataset *poDstDS;
    char **papszOptions = nullptr;
    poDstDS = poDriver->Create(outputFile, CELL_COUNT, CELL_COUNT, 3, GDT_UInt32, papszOptions);

    double adfGeoTransform[6] = { (double)bbox.minX, (double)xCellWidth, 0, (double)bbox.maxY, 0, -(double)yCellWidth };

    OGRSpatialReference oSRS;
    char *pszSRS_WKT = nullptr;
    GDALRasterBand *poBand;

    GUInt32 entryTimesRaster[CELL_COUNT * CELL_COUNT];
    GUInt32 exitTimesRaster[CELL_COUNT * CELL_COUNT];
    GUInt32 activeTimesRaster[CELL_COUNT * CELL_COUNT];

    poDstDS->SetGeoTransform( adfGeoTransform );

    std::cout << CELL_COUNT * CELL_COUNT << " " << timings.size() << "\n";

    for (int i = 0; i < CELL_COUNT * CELL_COUNT; i++){
        entryTimesRaster[i] = timings.at(i).firstTime;
        exitTimesRaster[i] = timings.at(i).lastTime;
        activeTimesRaster[i] = timings.at(i).lastTime - timings.at(i).firstTime;
    }

//    oSRS.SetWellKnownGeogCS("WGS84");
    oSRS.importFromEPSG(28992);
    oSRS.exportToWkt(&pszSRS_WKT);

    poDstDS->SetProjection(pszSRS_WKT);
    CPLFree(pszSRS_WKT);

    poBand = poDstDS->GetRasterBand(1);
    poBand->RasterIO(GF_Write, 0, 0, CELL_COUNT, CELL_COUNT, entryTimesRaster, CELL_COUNT, CELL_COUNT, GDT_UInt32, 0, 0);

    poBand = poDstDS->GetRasterBand(2);
    poBand->RasterIO(GF_Write, 0, 0, CELL_COUNT, CELL_COUNT, exitTimesRaster, CELL_COUNT, CELL_COUNT, GDT_UInt32, 0, 0);

    poBand = poDstDS->GetRasterBand(3);
    poBand->RasterIO(GF_Write, 0, 0, CELL_COUNT, CELL_COUNT, activeTimesRaster, CELL_COUNT, CELL_COUNT, GDT_UInt32, 0, 0);

    GDALClose((GDALDatasetH) poDstDS);


    lasreader->close();
    delete lasreader;

    return 0;
}
