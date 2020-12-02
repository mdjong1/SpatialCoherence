#include <iostream>
#include <utility>
#include <vector>
#include <cmath>
#include <string>
#include <filesystem>
#include <map>
#include <list>

#include "gdal_priv.h"
#include "cpl_conv.h"

#include "lasreader.hpp"

using namespace std;
namespace fs = filesystem;

struct Coordinate {
    double id, x, y;

    Coordinate(int paramID, double paramX, double paramY) : id(paramID), x(paramX), y(paramY) {}
};

struct Tile {
    int clusterId = -1;  // e.g. 37 (01 - 70)
    string rowInd;  // e.g. E (A - G)
    string rowSpec;  // e.g. N (always N or Z)
    int colId = -1;  // e.g. 1 (always 1 or 2)
    filesystem::path filepath; // path to file
};

string getTileNameFromTile(const Tile& inputTile) {
    return to_string(inputTile.clusterId) + inputTile.rowInd + inputTile.rowSpec + to_string(inputTile.colId);
}

string getNextTileName(Tile currentTile) {
    Tile nextTile;

    if (currentTile.clusterId + 1 <= 70) {

        if (currentTile.colId == 1) {
            nextTile.clusterId = currentTile.clusterId;
            nextTile.rowInd = currentTile.rowInd;
            nextTile.rowSpec = currentTile.rowSpec;
            nextTile.colId = 2;

        } else {
            nextTile.colId = 1;

            if (currentTile.rowSpec == "N") {
                nextTile.clusterId = currentTile.clusterId;
                nextTile.rowInd = currentTile.rowInd;
                nextTile.rowSpec = "Z";

            } else {
                nextTile.rowSpec = "N";

                if (currentTile.rowInd == "H") {  // Reached end of letter range, restarting
                    nextTile.rowInd = "A";
                    nextTile.clusterId = currentTile.clusterId + 1;

                } else {
                    char currentRowInd = currentTile.rowInd[0] + 1;  // ASCII value +1 == next letter in alphabet
                    nextTile.rowInd = currentRowInd;

                    nextTile.clusterId = currentTile.clusterId;
                }
            }
        }
    }

    return getTileNameFromTile(nextTile);
}

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
    if (returnValue % 10 == 9) { // Last digit is a 9
        returnValue++;
    }
    return returnValue;
}

Tile createTile(string inputData) {
    Tile outputTile;
    for (int i = inputData.size() - 1; i >= 0; i--) {
        if (isdigit(inputData[i]) && outputTile.colId == -1) {
            outputTile.colId = (int) inputData[i] - 48;  // Adjust for ASCII with -48
        } else if (!isdigit(inputData[i]) && outputTile.rowSpec.empty()) {
            outputTile.rowSpec = inputData[i];
        } else if (!isdigit(inputData[i]) && outputTile.rowInd == "") {
            outputTile.rowInd = inputData[i];
        } else if (isdigit(inputData[i]) && outputTile.clusterId == -1) {
            outputTile.clusterId = (int) inputData[i] - 48;
        } else if (isdigit(inputData[i]) && outputTile.clusterId != -1) {
            outputTile.clusterId = stoi(inputData.substr(i, 2));
        }
    }
    return outputTile;
}

string getTileNameFromPath(const filesystem::path &inputPath) {
    return inputPath.stem().generic_string().substr(inputPath.stem().generic_string().find('_') + 1, inputPath.stem().generic_string().size());
}

int main(int argc, char **argv) {
    if (argc != 7) {
        cout << "Invalid number of input arguments!\n";
        cout << "arg1: input folder, arg2: start tile, arg3: num tiles to process, arg4: output file, arg5: cell count, arg6: thinning_factor\n";
        return 0;
    }

    int streamTime = 0;

    const string inputFolder = argv[1];
    const string startTile = argv[2];
    const int numTilesToProcess = stoi(argv[3]);
    const char *outputFile = argv[4];
    const int CELL_COUNT = stoi(argv[5]); // stoi = cast to int
    const int THINNING_FACTOR = stoi(argv[6]);

    const int MAX_CELL_COUNT = CELL_COUNT + 2;

    GDALAllRegister();

    GUInt32 entryTimesRaster[MAX_CELL_COUNT * MAX_CELL_COUNT * numTilesToProcess];
    GUInt32 exitTimesRaster[MAX_CELL_COUNT * MAX_CELL_COUNT * numTilesToProcess];
    GUInt32 activeTimesRaster[MAX_CELL_COUNT * MAX_CELL_COUNT * numTilesToProcess];

    list<filesystem::path> eligibleFiles;

    for (const auto &entry : fs::directory_iterator(inputFolder)) {
        filesystem::path fileExtension = fs::path(entry.path()).extension();
        if (fileExtension == ".LAZ" || fileExtension == ".laz") {
            eligibleFiles.push_back(entry.path());
        }
    }

    map<string, Tile> tiles;

    for (const auto &eligibleFile : eligibleFiles) {
        Tile insertTile = createTile(eligibleFile.stem().generic_string());
        insertTile.filepath = eligibleFile;

        tiles.insert(pair<string, Tile>(getTileNameFromPath(eligibleFile), insertTile));
    }

    Tile currentTile = tiles[startTile];

    int cellWidth, cellHeight;
    int bboxMinX = INT_MAX;  // Ensure next value is always smaller
    int bboxMinY = INT_MAX;
    int bboxMaxX = INT_MIN;  // Ensure next value is always larger
    int bboxMaxY = INT_MIN;

    for (int tileNum = 0; tileNum < numTilesToProcess; tileNum++) {

        LASreadOpener lasreadopener;
        lasreadopener.set_file_name(currentTile.filepath.generic_string().c_str());  // filesystem::path -> string -> char*
        LASreader *lasreader = lasreadopener.open();

        // In AHN3 all corner points are integers, also round them up if precision is off (84999 -> 85000)
        const Bbox bbox = Bbox(roundUp(lasreader->get_min_x()), roundUp(lasreader->get_min_y()),
                               roundUp(lasreader->get_max_x()), roundUp(lasreader->get_max_y()));

        cout << "BBox boundaries: minX = " << bbox.minX << ", minY = " << bbox.minY << ", maxX = " << bbox.maxX
                  << ", maxY = " << bbox.maxY << "\n";

        if (bbox.minX < bboxMinX)
            bboxMinX = bbox.minX;
        if (bbox.minY < bboxMinY)
            bboxMinY = bbox.minY;
        if (bbox.maxX > bboxMaxX)
            bboxMaxX = bbox.maxX;
        if (bbox.maxY > bboxMaxY)
            bboxMaxY = bbox.maxY;

        const int numPoints = lasreader->npoints;

        cout << "Number of points: " << numPoints << "\n";

        const int xCellWidth = bbox.xDiff / CELL_COUNT;
        const int yCellWidth = bbox.yDiff / CELL_COUNT;

        cellWidth = xCellWidth;
        cellHeight = yCellWidth;

        Timings timings[MAX_CELL_COUNT][MAX_CELL_COUNT];

        int lastPercentage = -1;

        while (lasreader->read_point()) {

            streamTime++;

            if (streamTime % 20000 == 0) {
                int percentage = (int) round((double) streamTime / (double) numPoints * 100);
                if (percentage != lastPercentage) {
                    cout << percentage << "% done!\n";
                    lastPercentage = percentage;
                }
            }

            if (streamTime % THINNING_FACTOR == 0) {

                int xGridPos = MAX_CELL_COUNT - ((bbox.maxX - lasreader->point.get_x()) / xCellWidth);
                int yGridPos = (bbox.maxY - lasreader->point.get_y()) / yCellWidth;

                if (timings[xGridPos][yGridPos].firstTime == 0) {
                    timings[xGridPos][yGridPos].firstTime = streamTime;
                }

                timings[xGridPos][yGridPos].lastTime = streamTime;
            }
        }

        lasreader->close();
        delete lasreader;

        int rasterCellIndex = tileNum * MAX_CELL_COUNT;

        for (int y = 0; y < MAX_CELL_COUNT; y++) {
            for (int x = 0; x < MAX_CELL_COUNT; x++){
                entryTimesRaster[rasterCellIndex] = timings[x][y].firstTime;
                exitTimesRaster[rasterCellIndex] = timings[x][y].lastTime;
                activeTimesRaster[rasterCellIndex] = timings[x][y].lastTime - timings[x][y].firstTime;

                rasterCellIndex++;
            }
            rasterCellIndex += tileNum * MAX_CELL_COUNT;
        }

        currentTile = tiles[getNextTileName(currentTile)];
    }

    const int rasterXSize = (bboxMaxX - bboxMinX) / cellWidth;
    const int rasterYSize = (bboxMaxY - bboxMinY) / cellHeight;

    cout << "Writing GeoTIFF \n";

    const char *pszFormat = "GTiff";
    GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);

    if (poDriver == nullptr)
        exit(1);

    GDALDataset *poDstDS;
    char **papszOptions = nullptr;
    poDstDS = poDriver->Create(outputFile, rasterXSize, rasterYSize, 3, GDT_UInt32, papszOptions);

    double adfGeoTransform[6] = { (double)bboxMinX, (double)cellWidth, 0, (double)bboxMaxY, 0, -(double)cellHeight };

    OGRSpatialReference oSRS;
    char *pszSRS_WKT = nullptr;
    GDALRasterBand *poBand;

    poDstDS->SetGeoTransform( adfGeoTransform );

    //    oSRS.SetWellKnownGeogCS("WGS84");
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
