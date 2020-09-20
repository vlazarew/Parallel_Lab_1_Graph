#include "Portrait.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

const string portraitFilename = "D:/Универ/МГУ/5 курс/Параллельные вычисления/Лаба 1. Граф/GenerationOfPortraitOfMatrix/MakePortrait.txt";

Portrait createPortraitFromFile()
{
    ifstream input(portraitFilename);

    if (input.is_open()) {
        return readPortrainFromStream(input);
    } else {
        makeDefaultPortraitFile();
        createPortraitFromFile();
    }
}

Portrait readPortrainFromStream(ifstream& input)
{
    Portrait result;
    string line;
    int stringResults[4];

    try {
        while (getline(input, line)) {
            size_t pos = 0;
            string delimiter = " ";
            int index = 0;
            string token;

            while ((pos = line.find(delimiter)) != string::npos) {
                token = line.substr(0, pos);
                stringResults[index] = stoi(token);
                index++;
                line.erase(0, pos + delimiter.length());
            }
            stringResults[index] = stoi(line);
        }
    } catch (const std::exception&) {
    }

    result.setCountOfColumns(stringResults[0]);
    result.setCountOfRows(stringResults[1]);
    result.setK1(stringResults[2]);
    result.setK2(stringResults[3]);

    return result;
}

void makeDefaultPortraitFile()
{
    ofstream out;
    out.open(portraitFilename);

    if (out.is_open()) {
        out << "9 8 4 3" << endl;
    }

    cout << "Создан файл генерации портрета графа" << endl;
}