#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include <string>
#include <filesystem>
#include <map>
#include <list>
#include <algorithm>
#include <iterator>
#include <fstream>

#include "gdal_priv.h"
#include "cpl_conv.h"

#include "lasreader.hpp"

using namespace std;
namespace fs = filesystem;

struct Timings {
    int firstTime, lastTime;
};

int main(int argc, char **argv) {
    if (argc != 4) {
        cout << "Invalid number of input arguments!\n";
        cout << "arg1: input file, arg2: output file, arg3: cell count\n";
        return 0;
    }

    int streamTime = 0;

    const string inputFile = argv[1];
    const string outputFile = argv[2];
    const int CELL_COUNT = stoi(argv[3]); // stoi = cast to int

    const int MAX_CELL_COUNT = CELL_COUNT + 2;

    vector<string> inputFiles;

    ifstream inFile (inputFile);
    std::size_t lines_count =0;
    std::string line;
    while (std::getline(inFile , line)) {
        inputFiles.push_back(line);
        lines_count++;
    }

    const int numTilesToProcess = lines_count;

    GDALAllRegister();

    std::cout << "Files to process: " << numTilesToProcess << std::endl;

    GUInt32 entryTimesRaster[MAX_CELL_COUNT * MAX_CELL_COUNT * numTilesToProcess];
    GUInt32 exitTimesRaster[MAX_CELL_COUNT * MAX_CELL_COUNT * numTilesToProcess];
    GUInt32 activeTimesRaster[MAX_CELL_COUNT * MAX_CELL_COUNT * numTilesToProcess];

    int cellWidth, cellHeight;
    float bboxMinX = std::numeric_limits<float>::max();  // Ensure next value is always smaller
    float bboxMinY = std::numeric_limits<float>::max();
    float bboxMaxX = std::numeric_limits<float>::min();  // Ensure next value is always larger
    float bboxMaxY = std::numeric_limits<float>::min();

    Timings timings[MAX_CELL_COUNT][MAX_CELL_COUNT];

    for (const auto& currentFile: inputFiles) {

        LASreadOpener lasreadopener;
        lasreadopener.set_file_name(currentFile.c_str());  // filesystem::path -> string -> char*
        LASreader *lasreader = lasreadopener.open();

        float minX = lasreader->get_min_x();
        float maxX = lasreader->get_max_x();
        float minY = lasreader->get_min_y();
        float maxY = lasreader->get_max_y();

        if (minX < bboxMinX)
            bboxMinX = minX;
        if (minY < bboxMinY)
            bboxMinY = minY;
        if (maxX > bboxMaxX)
            bboxMaxX = maxX;
        if (maxY > bboxMaxY)
            bboxMaxY = maxY;

    }

    const int xCellWidth = (bboxMaxX - bboxMinX) / CELL_COUNT;
    const int yCellWidth = (bboxMaxY - bboxMinY) / CELL_COUNT;

    const int rasterXSize = (bboxMaxX - bboxMinX) / xCellWidth;
    const int rasterYSize = (bboxMaxY - bboxMinY) / yCellWidth;

    std::cout << "XSize: " << rasterXSize << " | YSize: " << rasterYSize << std::endl;

    for (const auto& currentFile: inputFiles) {

        std::cout << "Current file = " << currentFile << std::endl;

        LASreadOpener lasreadopener;
        lasreadopener.set_file_name(currentFile.c_str());  // filesystem::path -> string -> char*
        LASreader *lasreader = lasreadopener.open();

        const int numPoints = lasreader->npoints;

        cout << "Number of points: " << numPoints << "\n";

        int lastPercentage = -1;
        int pointCount = 0;

        while (lasreader->read_point()) {

            streamTime++;
            pointCount++;

            if (pointCount % 20000 == 0) {
                int percentage = (int) round(((double) pointCount / (double) numPoints) * 100);
                if (percentage != lastPercentage) {
                    cout << percentage << "% done!\n";
                    lastPercentage = percentage;
                }
            }

            int xGridPos = MAX_CELL_COUNT - ((bboxMaxX - lasreader->point.get_x()) / xCellWidth);
            int yGridPos = (bboxMaxY - lasreader->point.get_y()) / yCellWidth;

            if (timings[xGridPos][yGridPos].firstTime == 0) {
                timings[xGridPos][yGridPos].firstTime = streamTime;
            }

            timings[xGridPos][yGridPos].lastTime = streamTime;
        }

        lasreader->close();
        delete lasreader;
    }

    int i = 0;
    for (int y = 0; y < CELL_COUNT; y++) {
        for (int x = 0; x < CELL_COUNT; x++){
            entryTimesRaster[i] = timings[x][y].firstTime;
            exitTimesRaster[i] = timings[x][y].lastTime;
            activeTimesRaster[i] = timings[x][y].lastTime - timings[x][y].firstTime;

            i++;
        }
    }

    cout << "Writing GeoTIFF \n";

    const char *pszFormat = "GTiff";
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

    if (poDriver == nullptr)
        exit(1);

    GDALDataset *poDstDS;
    char **papszOptions = nullptr;
    poDstDS = poDriver->Create(outputFile.c_str(), rasterXSize, rasterYSize, 3, GDT_UInt32, papszOptions);

    double adfGeoTransform[6] = { (double)bboxMinX, (double)xCellWidth, 0, (double)bboxMaxY, 0, -(double)yCellWidth };

    OGRSpatialReference oSRS;
    char *pszSRS_WKT = nullptr;
    GDALRasterBand *poBand;

    poDstDS->SetGeoTransform( adfGeoTransform );

    oSRS.importFromEPSG(28992);
    oSRS.exportToWkt(&pszSRS_WKT);

    poDstDS->SetProjection(pszSRS_WKT);
    CPLFree(pszSRS_WKT);

    poBand = poDstDS->GetRasterBand(1);
    poBand->RasterIO(GF_Write, 0, 0, rasterXSize, rasterYSize, entryTimesRaster, rasterXSize, rasterYSize, GDT_UInt32, 0, 0);

    poBand = poDstDS->GetRasterBand(2);
    poBand->RasterIO(GF_Write, 0, 0, rasterXSize, rasterYSize, exitTimesRaster, rasterXSize, rasterYSize, GDT_UInt32, 0, 0);

    poBand = poDstDS->GetRasterBand(3);
    poBand->RasterIO(GF_Write, 0, 0, rasterXSize, rasterYSize, activeTimesRaster, rasterXSize, rasterYSize, GDT_UInt32, 0, 0);

    GDALClose((GDALDatasetH) poDstDS);

    return 0;
}
