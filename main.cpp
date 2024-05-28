#include <iostream>
#include <vector>
#include <string>
#include <iomanip> 

using namespace std;

int alpha(char a, char b) {
  if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A') ||
      (a == 'C' && b == 'G') || (a == 'G' && b == 'C')) {
    return -1; 
  }
  return 0; 
}

int alpha_2(char a, char b) {
  if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A')){
    return -4;
  }
  else if ((a == 'C' && b == 'G') || (a == 'G' && b == 'C')) {
    return -5; 
  }
  return 0; 
}


int sec_ARN(const string& secuencia, vector<vector<int>>& E, vector<vector<int>>& P) {
  int L = secuencia.length();

  E.assign(L, vector<int>(L, 0));
  P.assign(L, vector<int>(L, -1));

  for (int d = 2; d <= L; ++d) {
    for (int i = 0; i <= L - d; ++i) {
      int j = i + d - 1;

      int costo = E[i+1][j-1] + alpha(secuencia[i], secuencia[j]);
      if (costo < E[i][j-1]) {
        E[i][j] = costo;
        P[i][j] = -2;
      } else {
        E[i][j] = E[i][j-1];
        P[i][j] = -1;
      }

      for (int k = i + 1; k < j; ++k) {
        if (E[i][j] > E[i][k] + E[k+1][j]) {
          E[i][j] = E[i][k] + E[k+1][j];
          P[i][j] = k;
        }
      }
    }
  }

  return E[0][L-1];
}

void imprimirMatriz(const vector<vector<int>>& E) {
  int L = E.size();
  for (int i = 0; i < L; ++i) {
    for (int j = 0; j < L; ++j) {
      cout << setw(3) << E[i][j] << " ";
    }
    cout << endl;
  }
}

int main() {
  string secuencia = "GGAAAUCC";
  // string secuencia = "ACUCGAUUCCGAG";
  vector<vector<int>> E, P;
  int resultado = sec_ARN(secuencia, E, P);
  cout << "Matriz de scores:" << endl;
  imprimirMatriz(E);
  cout << "Costo mínimo de emparejamiento: " << resultado << endl;
  return 0;
}