#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>

#include "gdal_priv.h"
#include "cpl_conv.h"

#include "lasreader.hpp"

const int cellCount = 128;

struct Coordinate {
    double id, x, y;

    Coordinate(int paramID, double paramX, double paramY) : id(paramID), x(paramX), y(paramY) {}
};

struct Timings {
    int firstTime, lastTime;

    Timings(int first, int last) : firstTime(first), lastTime(last) {}
};

int main(int argc, char **argv) {
    if (argc <= 1) {
        std::cout << "Not enough input arguments, arg1 must be file name";
        return 0;
    } else if (argc > 2) {
        std::cout << "Too many input arguments!";
        return 0;
    }

    const char *inputFile = argv[1];

    LASreadOpener lasreadopener;
    lasreadopener.set_file_name(inputFile);
    LASreader *lasreader = lasreadopener.open();

    // In AHN3 all corner points are integers
    const int bbox[4] = {(int) lasreader->get_min_x(), (int) lasreader->get_min_y(),
                         (int) lasreader->get_max_x(), (int) lasreader->get_max_y()};

    std::cout << "BBox boundaries: minX = " << bbox[0] << ", minY = " << bbox[1] << ", maxX = " << bbox[2]
              << ", maxY = " << bbox[3] << "\n";

    const int numPoints = lasreader->npoints;

    std::cout << "Amount of points: " << numPoints << "\n";

    const int xDiff = bbox[2] - bbox[0];
    const int yDiff = bbox[3] - bbox[1];

    const int xCellWidth = xDiff / cellCount;
    const int yCellWidth = yDiff / cellCount;

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

    while (lasreader->read_point()) {

        if (count % 20000 == 0) {
            double percentage = round((double) count / (double) numPoints * 100);
            std::cout << percentage << "% done!\n";
        }

        for (Coordinate corner : leftTopCorners) {

            if (lasreader->point.inside_rectangle(corner.x, corner.y, corner.x + xCellWidth, corner.y + yCellWidth)) {

//                std::cout << lasreader->point.get_x() << ", " << lasreader->point.get_y() << " inside " << corner.x << ", " << corner.y << ", " << corner.x + xCellWidth << ", " << corner.y + yCellWidth << "\n";
//                std::cout << "corner id = " << corner.id << "\n";

                int currentEpoch = std::chrono::duration_cast<std::chrono::seconds>(
                        std::chrono::system_clock::now().time_since_epoch()).count();

                try { // Vector exists, update last entry time
                    timings.at(corner.id).lastTime = currentEpoch;
                }
                catch (...) { // Vector doesn't exist, set first + last time and push back
                    timings.emplace_back(currentEpoch, currentEpoch);
                }

                break;
            }
        }

        count++;

//        if (count > 50000){
//            break;
//        }
    }

//    for (Timings result : timings) {
////        std::cout << result.firstTime << " -> " << result.lastTime << "\n";
//        std::cout << result.lastTime - result.firstTime << "\n";
//    }


    GDALAllRegister();

    const char *pszDstFilename = "test4.tif";

    const char *pszFormat = "GTiff";
    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

    if (poDriver == nullptr)
        exit(1);

    GDALDataset *poDstDS;
    char **papszOptions = nullptr;
    poDstDS = poDriver->Create(pszDstFilename, cellCount, cellCount, 1, GDT_Byte, papszOptions);

    double adfGeoTransform[6] = { (double)bbox[0], (double)xCellWidth, 0, (double)bbox[3], 0, -(double)yCellWidth };

    OGRSpatialReference oSRS;
    char *pszSRS_WKT = nullptr;
    GDALRasterBand *poBand;

    GByte abyRaster[cellCount * cellCount];

    poDstDS->SetGeoTransform( adfGeoTransform );

    std::cout << cellCount * cellCount << " " << timings.size();

    for (int i = 0; i < cellCount * cellCount; i++){
        abyRaster[i] = timings.at(i).lastTime - timings.at(i).firstTime;
    }

//    oSRS.SetWellKnownGeogCS("WGS84");
    oSRS.importFromEPSG(28992);
    oSRS.exportToWkt(&pszSRS_WKT);

    poDstDS->SetProjection(pszSRS_WKT);
    CPLFree(pszSRS_WKT);

    poBand = poDstDS->GetRasterBand(1);

    poBand->RasterIO(GF_Write, 0, 0, cellCount, cellCount, abyRaster, cellCount, cellCount, GDT_Byte, 0, 0);

    GDALClose((GDALDatasetH) poDstDS);


    lasreader->close();
    delete lasreader;

    return 0;
}
