#pragma once

#include "const.cuh"

int CHECK(int n);

template<typename T>
void SAY_LOG(T s);

template<typename T>
void stp(T str);

template<typename T>
string vectorToStr(vector<T> v, string sp = ";");

template<typename T>
void printVector(ostream &output, vector<T> v, string sp1="(", string sp2=";", string sp3=")\n");

template<typename T>
T sumVector(vector<T> v);

template <typename T>
string toString(T val);

template<typename T>
T fromString(const string& s);

string toLower(string s);
string toUpper(string s);
void printX3(ostream &output, double3 r, string sp1="(", string sp2=";", string sp3=")\n");
void printX3(ostream &output, int3 r, string sp1="(", string sp2=";", string sp3=")\n");
string toStr(double3 d, string sp1="(", string sp2=";", string sp3=")");
string toStr(int3 d, string sp1="(", string sp2=";", string sp3=")");
vector<double> d3ToV(double3 v);

void time_progress(time_t real_start_t, time_t curr_t, double done_part, string proc_name);
