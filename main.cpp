#include <iostream>
#include <vector>
#include <chrono>
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

int nearestThousandth(F64 number) {
    return (int) round(number / 1000) * 1000;
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

    // In AHN3 all corner points are integers, round them to nearest thousandth (84999 -> 85000)
    const Bbox bbox = Bbox(nearestThousandth(lasreader->get_min_x()), nearestThousandth(lasreader->get_min_y()),
                           nearestThousandth(lasreader->get_max_x()), nearestThousandth(lasreader->get_max_y()));

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

    count = 0;
    const int startEpoch = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    while (lasreader->read_point()) {

        if (count % 20000 == 0) {
            double percentage = round((double) count / (double) numPoints * 100);
            std::cout << percentage << "% done!\n";
        }

        if (count % THINNING_FACTOR == 0) {

            for (Coordinate corner : leftTopCorners) {

                if (lasreader->point.inside_rectangle(corner.x, corner.y, corner.x + xCellWidth, corner.y + yCellWidth)) {

//                    std::cout << lasreader->point.get_x() << ", " << lasreader->point.get_y() << " inside " << corner.x << ", " << corner.y << ", " << corner.x + xCellWidth << ", " << corner.y + yCellWidth << "\n";
//                    std::cout << "corner id = " << corner.id << "\n";

                    int currentEpoch = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();

                    try { // Vector exists, update last entry time
                        timings.at(corner.id).lastTime = currentEpoch;
                    }
                    catch (...) { // Vector doesn't exist, set first + last time and push back
                        timings.emplace_back(currentEpoch, currentEpoch);
                    }

                    break;
                }
            }
        }

        count++;

    }

    std::cout << "Writing GeoTIFF \n";

    GDALAllRegister();

    const char *pszFormat = "GTiff";
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

    if (poDriver == nullptr)
        exit(1);

    GDALDataset *poDstDS;
    char **papszOptions = nullptr;
    poDstDS = poDriver->Create(outputFile, CELL_COUNT, CELL_COUNT, 3, GDT_Byte, papszOptions);

    double adfGeoTransform[6] = { (double)bbox.minX, (double)xCellWidth, 0, (double)bbox.maxY, 0, -(double)yCellWidth };

    OGRSpatialReference oSRS;
    char *pszSRS_WKT = nullptr;
    GDALRasterBand *poBand;

    GByte entryTimesRaster[CELL_COUNT * CELL_COUNT];
    GByte exitTimesRaster[CELL_COUNT * CELL_COUNT];
    GByte activeTimesRaster[CELL_COUNT * CELL_COUNT];

    poDstDS->SetGeoTransform( adfGeoTransform );

    std::cout << CELL_COUNT * CELL_COUNT << " " << timings.size() << "\n";

    for (int i = 0; i < CELL_COUNT * CELL_COUNT; i++){
        entryTimesRaster[i] = timings.at(i).firstTime - startEpoch;
        exitTimesRaster[i] = timings.at(i).lastTime - startEpoch;
        activeTimesRaster[i] = timings.at(i).lastTime - timings.at(i).firstTime;
    }

//    oSRS.SetWellKnownGeogCS("WGS84");
    oSRS.importFromEPSG(28992);
    oSRS.exportToWkt(&pszSRS_WKT);

    poDstDS->SetProjection(pszSRS_WKT);
    CPLFree(pszSRS_WKT);

    poBand = poDstDS->GetRasterBand(1);
    poBand->RasterIO(GF_Write, 0, 0, CELL_COUNT, CELL_COUNT, entryTimesRaster, CELL_COUNT, CELL_COUNT, GDT_Byte, 0, 0);

    poBand = poDstDS->GetRasterBand(2);
    poBand->RasterIO(GF_Write, 0, 0, CELL_COUNT, CELL_COUNT, exitTimesRaster, CELL_COUNT, CELL_COUNT, GDT_Byte, 0, 0);

    poBand = poDstDS->GetRasterBand(3);
    poBand->RasterIO(GF_Write, 0, 0, CELL_COUNT, CELL_COUNT, activeTimesRaster, CELL_COUNT, CELL_COUNT, GDT_Byte, 0, 0);

    GDALClose((GDALDatasetH) poDstDS);


    lasreader->close();
    delete lasreader;

    return 0;
}
