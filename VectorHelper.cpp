#include "VectorHelper.h"


std::vector<uint8_t> getRandomVector(int size) {
    std::vector<uint8_t> ret(size);

    for (int i = 0; i < size; i++) {
        ret[i] = rand() % 2;
    }

    return ret;
}

std::string getStringHexFromVector(std::vector<uint8_t> data) {
    std::string ret = "";
    for (int i = data.size() / 4 - 1; i >= 0; i--) {
        uint8_t hex = 0;
        for (int j = 0; j < 4; j++) {
            hex = (hex << 1) | data[i * 4 + 3 - j];
        }

        if (hex < 10)
            ret.push_back(hex + '0');
        else
            ret.push_back(hex + 'A' - 10);
    }

    return ret;
}

std::string getStringBinFromVector(std::vector<uint8_t> data) {
    std::string ret = "";
    for (int i = 0; i < data.size(); i++) {
    //for (int i = data.size() - 1; i >= 0; i--) {
        ret.push_back(data[i] + '0');
        //ret.push_back(',');
    }

    return ret;
}


std::vector<uint8_t> getVectorFromStringHex(std::string s) {
    std::vector<uint8_t> ret;

    //for (int i = 0; i < s.length(); i++) {
    for (int i = s.length() - 1; i >= 0; i--) {
        uint8_t halfByte;

        if (s[i] >= 'a' && s[i] <= 'f') {
            halfByte = 10 + s[i] - 'a';
        }
        else
            if (s[i] >= 'A' && s[i] <= 'F') {
                halfByte = 10 + s[i] - 'A';
            }
            else
                if (s[i] >= '0' && s[i] <= '9') {
                    halfByte = s[i] - '0';
                }
                else
                    continue;

        for (int j = 0; j < 4; j++) {
            ret.push_back((halfByte >> j) & 1);
        }
    }

    return ret;
}

std::vector<uint8_t> getVectorFromStringBin(std::string s) {
    std::vector<uint8_t> ret;

    //for (int i = 0; i < s.length(); i++) {
    for (int i = s.length() - 1; i >= 0; i--) {
        if (s[i] == '0')
            ret.push_back(0);
        else
            if (s[i] == '1')
                ret.push_back(1);
    }

    return ret;
}

std::vector<uint8_t> getVectorFromArray(uint8_t* array, int size) {
    std::vector<uint8_t> ret(size);

    for (int i = 0; i < size; i++) {
        ret[i] = array[i];
    }

    return ret;
}

int getErrors(std::vector<uint8_t> a, std::vector<uint8_t> b) {
    int err = 0;
    for (int i=0; i< std::min(a.size(), b.size()); i++) {
        if (a[i] != b[i])
            err++;
    }
    return  err;
}
