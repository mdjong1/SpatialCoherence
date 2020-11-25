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
    int CELL_COUNT = std::stoi(argv[3]); // std::stoi = cast to int
    const int THINNING_FACTOR = std::stoi(argv[4]);

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

    CELL_COUNT += 2;

    Timings timings[CELL_COUNT][CELL_COUNT];

    int streamTime = 0;
    int lastPercentage = -1;

    while (lasreader->read_point()) {

        streamTime++;

        if (streamTime % 20000 == 0) {
            int percentage = (int) round((double) streamTime / (double) numPoints * 100);
            if (percentage != lastPercentage) {
                std::cout << percentage << "% done!\n";
                lastPercentage = percentage;
            }
        }

        if (streamTime % THINNING_FACTOR == 0) {

            int xGridPos = CELL_COUNT - ((bbox.maxX - lasreader->point.get_x()) / xCellWidth);
            int yGridPos = (bbox.maxY - lasreader->point.get_y()) / yCellWidth;

            if (timings[xGridPos][yGridPos].firstTime == 0) {
                timings[xGridPos][yGridPos].firstTime = streamTime;
            }

            timings[xGridPos][yGridPos].lastTime = streamTime;

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

    int i = 0;
    for (int x = 0; x < CELL_COUNT; x++){
        for (int y = 0; y < CELL_COUNT; y++) {
            entryTimesRaster[i] = timings[x][y].firstTime;
            exitTimesRaster[i] = timings[x][y].lastTime;
            activeTimesRaster[i] = timings[x][y].lastTime - timings[x][y].firstTime;

            i++;
        }
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
