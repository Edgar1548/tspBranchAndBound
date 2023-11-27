#include <bits/stdc++.h>
#include <fstream>
#include <sstream>
#include <omp.h>


using namespace std;

// Función para calcular la distancia
double calculateDistance(int x1, int y1, int x2, int y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

// Función para convertir coordenadas en una matriz de distancias
vector<vector<double>> convertToMatrix(const vector<pair<int, int>>& coordinates) {
    int n = coordinates.size();
    vector<vector<double>> distanceMatrix(n, vector<double>(n));

    #pragma omp parallel for collapse(2)
    for (int a = 0; a < n; ++a) {
        for (int b = 0; b < n; ++b) {
            if (a == b) distanceMatrix[a][b] = INT_MAX;
            else distanceMatrix[a][b] = calculateDistance(coordinates[a].first, coordinates[a].second, coordinates[b].first, coordinates[b].second);
        }
    }
    return distanceMatrix;
}


// Variables globales ajustadas para trabajar con matrices dinámicas
vector<vector<double>> adj;
vector<int> final_path;
vector<bool> visited;
double final_res = DBL_MAX;
int N;

void copyToFinal(const vector<int>& curr_path) {
    for (int i = 0; i < N; i++)
        final_path[i] = curr_path[i];
    final_path[N] = curr_path[0];
}

double firstMin(int i) {
    double min = DBL_MAX;
    for (int k = 0; k < N; k++)
        if (adj[i][k] < min && i != k)
            min = adj[i][k];
    return min;
}


double secondMin(int i) {
    double first = DBL_MAX, second = DBL_MAX;
    for (int j = 0; j < N; j++) {
        if (i == j)
            continue;

        if (adj[i][j] <= first) {
            second = first;
            first = adj[i][j];
        } else if (adj[i][j] <= second && adj[i][j] != first)
            second = adj[i][j];
    }
    return second;
}




void TSPRec(double curr_bound, double curr_weight, int level, vector<int>& curr_path, vector<bool>& visited) {
    if (level == N) {
        if (adj[curr_path[level - 1]][curr_path[0]] != INT_MAX) {
            double curr_res = curr_weight + adj[curr_path[level - 1]][curr_path[0]];
            #pragma omp critical
            {
                if (curr_res < final_res) {
                    copyToFinal(curr_path);
                    final_res = curr_res;
                }
            }
        }
        return;
    }

    for (int i = 0; i < N; i++) {
        if (adj[curr_path[level - 1]][i] != INT_MAX && !visited[i]) {
            double temp_curr_weight = curr_weight + adj[curr_path[level - 1]][i];
            double temp_curr_bound = curr_bound;

            if (level == 1)
                temp_curr_bound -= ((firstMin(curr_path[level - 1]) + firstMin(i)) / 2);
            else
                temp_curr_bound -= ((secondMin(curr_path[level - 1]) + firstMin(i)) / 2);

            if (temp_curr_bound + temp_curr_weight < final_res) {
                visited[i] = true;
                curr_path[level] = i;

                TSPRec(temp_curr_bound, temp_curr_weight, level + 1, curr_path, visited);

                // Restaurar estado para la siguiente iteración
                visited[i] = false;
            }
        }
    }
}






void TSP() {
    vector<int> curr_path(N + 1, -1); // Almacena el camino actual en la búsqueda
    double curr_bound = 0;            // Almacena el límite inferior del camino
    std::fill(visited.begin(), visited.end(), false); // Inicializa todos los nodos como no visitados

    // Calcula el límite inferior inicial para el camino
    for (int i = 0; i < N; i++)
        curr_bound += (firstMin(i) + secondMin(i));

    // Redondea el límite inferior al número entero más cercano
    curr_bound = fmod(curr_bound, 2.0) > 0.0 ? curr_bound / 2 + 1 : curr_bound / 2;

    // Comienza desde el primer nodo
    visited[0] = true;
    curr_path[0] = 0;

    // Llama a la función recursiva TSPRec para construir el camino
    TSPRec(curr_bound, 0, 1, curr_path, visited);
}



// Funciones copyToFinal, firstMin, secondMin, TSPRec, TSP...
// Estas funciones deben modificarse para trabajar con 'vector<vector<double>> adj' y 'vector<int> final_path', etc.

// Función principal
int main(int argc, char *argv[]) {
    omp_set_num_threads(12);
    vector<pair<int, int>> coords;
    if (argc != 2) {
        cout << "No argument provided." << endl;
        return -1;
    } else {
        string line;
        ifstream myfile(argv[1]);
        if (!myfile) {
            cout << "Cannot open the file." << endl;
            return -1;
        }

        while (getline(myfile, line)) {
            istringstream iss(line);
            int nodo, x, y;
            if (!(iss >> nodo >> x >> y)) {
                continue; // Manejar errores de formato
            }
            coords.push_back({x, y});
        }
        myfile.close();
    }

    // Iniciar el temporizador
    auto start = chrono::high_resolution_clock::now();

    
    adj = convertToMatrix(coords);
    N = adj.size();
    final_path.resize(N + 1);
    visited.resize(N, false);



    // Llamada a la función TSP
    TSP();

    // Detener el temporizador
    auto stop = chrono::high_resolution_clock::now();

    // Calcular la duración
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);

    cout << "Minimum cost : " << final_res << endl;
    cout << "Path Taken : ";
    for (int i = 0; i <= N; i++)
        cout << final_path[i] << " ";
    cout << endl;

    // Mostrar la duración
    cout << "Time taken by function: " << duration.count() << " milliseconds" << endl;

    return 0;
}
