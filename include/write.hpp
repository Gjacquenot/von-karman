#ifndef WRITE_HPP
#define WRITE_HPP

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "misc.hpp"

using namespace std;

template <typename... Args>

// @brief Write the arguments to a temporary file
// @param args arguments to be written to the file
void write_tmp_file(Args... args) {
  ofstream file;
  file.open("data/tmp_animation.txt");
  // write all the arguments to the file one by one (separated by a space)
  ((file << args << " "), ...);
  file.close();
}

// @brief Print a string and a time in a human-readable format (in the form: "string xx.xx units")
// @param str string to be printed
// @param time time to be printed
void print(string str, int64_t time) {
  // default time is in microseconds
  string microsec = "\u03BCs";
  string millisec = "ms";
  string sec = "s";
  string min = "min";
  string hour = "h";
  double time_d;
  if (time > 10000) {
    time /= 1000;  // now time is in milliseconds
    if (time > 10000) {
      time /= 1000;  // now time is in seconds
      if (time > 100) {
        time_d = (double)time;
        time_d /= 60;  // now time is in minutes
        if (time_d > 100) {
          time_d /= 60;  // now time is in hours
          cout << str << time_d << " " << hour << endl;
          return;
        }
        cout << str << time_d << " " << min << endl;
        return;
      }
      cout << str << time << " " << sec << endl;
      return;
    }
    cout << str << time << " " << millisec << endl;
    return;
  }
  cout << str << time << " " << microsec << endl;
  return;
}

// @brief Write the euclidean norm of the velocity field to a file
// @param file file to write the solution
// @param u velocity in x direction
// @param v velocity in y direction
// @param t time at which the solution is written
// @param prm parameters of the simulation (dx, dy, dt, etc.)
void write_sol(ofstream &file, double *u, double *v, double t, Prm prm) {
  file << t << endl;
  for (int i = 1; i < prm.NX - 1; i++) {
    for (int j = 1; j < prm.NY - 1; j++) {
      file << sqrt(U(i, j) * U(i, j) + V(i, j) * V(i, j)) << " ";
    }
    file << endl;
  }
  file << endl;
}

// @brief Write the vorticity field to a file
// @param file file to write the solution
// @param w vorticity
// @param t time at which the solution is written
// @param prm parameters of the simulation (dx, dy, dt, etc.)
void write_sol_w(ofstream &file, double *w, double t, Prm prm) {
  file << t << endl;
  for (int i = 1; i < prm.NX - 1; i++) {
    for (int j = 1; j < prm.NY - 1; j++) {
      file << W(i, j) << " ";
    }
    file << endl;
  }
  file << endl;
}
#endif  // WRITE_HPP