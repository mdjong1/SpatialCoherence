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

    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(inputFile);
    LASreader *lasreader = lasreadopener.open();

    // In AHN3 all corner points are integers
    const int bbox[4] = {(int) lasreader->get_min_x(), (int) lasreader->get_min_y(),
                         (int) lasreader->get_max_x(), (int) lasreader->get_max_y()};

    std::cout << "BBox boundaries: minX = " << bbox[0] << ", minY = " << bbox[1] << ", maxX = " << bbox[2]
              << ", maxY = " << bbox[3] << "\n";

    const int numPoints = lasreader->npoints;

    std::cout << "Number of points: " << numPoints << "\n";

    const int xDiff = bbox[2] - bbox[0];
    const int yDiff = bbox[3] - bbox[1];

    const int xCellWidth = xDiff / CELL_COUNT;
    const int yCellWidth = yDiff / CELL_COUNT;

    std::vector<Coordinate> leftTopCorners;

    int count = 0;

    for (int x = bbox[0]; x < bbox[2]; x += xCellWidth) {
        for (int y = bbox[1]; y < bbox[3]; y += yCellWidth) {
            leftTopCorners.emplace_back(count, x, y);
            count++;
        }
    }

    std::vector<Timings> timings;

    count = 0;
    const int startEpoch = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    while (lasreader->read_point()) {

        if (count % 20000 == 0) {
            double percentage = round((double) count / (double) numPoints * 100);
            std::cout << percentage << "% done!\n";
        }

        if (count % THINNING_FACTOR == 0) {

            for (Coordinate corner : leftTopCorners) {

                if (lasreader->point.inside_rectangle(corner.x, corner.y, corner.x + xCellWidth,
                                                      corner.y + yCellWidth)) {

//                std::cout << lasreader->point.get_x() << ", " << lasreader->point.get_y() << " inside " << corner.x << ", " << corner.y << ", " << corner.x + xCellWidth << ", " << corner.y + yCellWidth << "\n";
//                std::cout << "corner id = " << corner.id << "\n";

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
    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

    if (poDriver == nullptr)
        exit(1);

    GDALDataset *poDstDS;
    char **papszOptions = nullptr;
    poDstDS = poDriver->Create(outputFile, CELL_COUNT, CELL_COUNT, 2, GDT_Byte, papszOptions);

    double adfGeoTransform[6] = { (double)bbox[0], (double)xCellWidth, 0, (double)bbox[3], 0, -(double)yCellWidth };

    OGRSpatialReference oSRS;
    char *pszSRS_WKT = nullptr;
    GDALRasterBand *poBand;

    GByte entryTimesRaster[CELL_COUNT * CELL_COUNT];
    GByte exitTimesRaster[CELL_COUNT * CELL_COUNT];

    poDstDS->SetGeoTransform( adfGeoTransform );

    std::cout << CELL_COUNT * CELL_COUNT << " " << timings.size() << "\n";

    for (int i = 0; i < CELL_COUNT * CELL_COUNT; i++){
        entryTimesRaster[i] = timings.at(i).firstTime - startEpoch;
        exitTimesRaster[i] = timings.at(i).lastTime - startEpoch;
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

    GDALClose((GDALDatasetH) poDstDS);


    lasreader->close();
    delete lasreader;

    return 0;
}
