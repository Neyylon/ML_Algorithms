
#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>

class Point2D {
private:
    std::pair<double, double> data;
public:
    Point2D() = default;
    Point2D(double Data1stPart, double Data2ndPart) noexcept {
        data.first = Data1stPart;
        data.second = Data2ndPart;
    }
    Point2D(std::pair<double, double> RecieveData) noexcept {
        data = RecieveData;
    }
    Point2D(std::initializer_list<double> Data) {
        if (Data.size() != 2)
            throw std::invalid_argument("Invalid argument. Spacify only 2 numbers.");
        auto it = Data.begin();
        data.first = *it++;
        data.second = *it;
    }

    inline constexpr double X()const noexcept { return data.first; }
    inline constexpr double Y()const noexcept { return data.second; }
}; // SML

class PointND {
private:
    size_t numberOfElements;
    std::vector<double> element;
public:
    PointND() = delete;
    PointND(std::initializer_list<double> coordinates) {
        for (auto& el : coordinates) {
            element.push_back(el);
        }
        numberOfElements = element.size();
    }

    // inline constexpr double fir()const noexcept { return ; }
};

std::ostream& operator<<(std::ostream& stream, const Point2D& point) {
    stream << point.X() << ' ' << point.Y();
    return stream;
}

std::ostream& operator<<(std::ostream& stream, const std::vector<Point2D>& points) {
    for (const auto& point : points) {
        stream << point << std::endl;
    }
    return stream;
}
//constexpr calculating during compilation time
//inline make the function faster

void Point2DUnitTest();
std::vector<Point2D> ReadFile(std::string filename, char delimeter, size_t X_col_index = 0, size_t Y_col_index = 1);

std::pair<double, double> FindCoeficient(std::vector<Point2D> points);

class Matrix {
private:
    std::vector<double> values;
    size_t rows;
    size_t columns;
public:
    Matrix() = delete;
    Matrix( size_t rows, size_t columns) {
        this->rows = rows;
        this->columns = columns;
    }
    //Matrix(std::string filename, char delimeter, );
    Matrix(std::initializer_list < std::vector<double>> matrix) {
        for (auto& row : matrix) 
            for (auto item : row)
                values.push_back(item);       
        rows = matrix.size();
        columns = matrix.begin()->size();
    }


    constexpr bool IsMultipliable(const Matrix& other)const noexcept { return columns == other.rows; }

    double& at(size_t row, size_t column) {
        if (row >= rows || column >= columns)
            throw std::invalid_argument("Choose row and column which less or equal to rows and columns of the matrix.");
        return values[row * columns + column];
    }

    double at(size_t row, size_t column)const {
        if (row >= rows || column >= columns)
            throw std::invalid_argument("Choose row and column which less or equal to rows and columns of the matrix.");
        return values[row * columns + column];
    }

    Matrix& Transpose() {
        for (size_t i{ 1 }; i < values.size(); ++i) {
            for (size_t j{ 0 }; j < i; ++j) {
                values[i * columns + j] ^= values[j * columns + i];
                values[j * columns + i] ^= values[i * columns + j];
                values[i * columns + j] ^= values[j * columns + i];
            }
        }
        return *this;
    }

    Matrix& operator*(const double number) {
        std::for_each(values.begin(), values.end(), [number](double& elem) {elem *= number; });
        return *this;
    }

    Matrix& operator*(const Matrix& other) {
        if (!IsMultipliable(other))
            throw std::invalid_argument("Write correct matrix");

            std::vector<double> Final;          
            for (size_t i{ 0 }; i < rows; ++i) {
                for (size_t j{ 0 }; j < other.columns; ++j) {
                    double accumulate{ 0.0 };
                    for (size_t k{ 0 }; k < other.rows; ++k) {
                        accumulate += values[i * columns + k] * other.values[k * other.columns + j];
                    }
                    Final.push_back(accumulate);
                }
            }

            values = std::move(Final);
            columns = other.columns;
            return *this;
    }

    double Determinant() const{
        if (rows != columns)
            throw std::invalid_argument("Set a square matrix.");
        if (values.size() == 1)
            return values[0];
        if (values.size() == 4)
            return values[0] * values[3] - values[1] * values[2];

        double accumulate{0.0};
        for (size_t i{ 0 }; i < rows; ++i) { 
            accumulate += Minor(i, 0).Determinant() * values[i*columns] * static_cast<int>(pow(-1, i + 2));
            
        }
        
        return accumulate;
    }

    Matrix Minor(size_t i, size_t j) const{
        Matrix temp(rows-1, columns-1);
        for (size_t k_of_i{0}; k_of_i < rows; ++k_of_i) {
            for (size_t t_of_j{0}; t_of_j < columns; ++t_of_j) {
                if(k_of_i != i && t_of_j != j)
                    temp.values.push_back(values[k_of_i* columns + t_of_j]);
            }
        }
        return temp;
    }

    Matrix Invertion() {
        Matrix temp(rows, columns);
        double det = this->Determinant();
        for (size_t i{ 0 }; i < rows; ++i) {
            for (size_t j{ 0 }; j < columns; ++j) {
                temp.values.push_back(this->Minor(i, j).Determinant() * pow(-1, i + j) / det);
            }
        }
        return temp.Transpose();
    }

  /*  Matrix E_matrix() {
        Matrix temp((rows * columns), 0);
        size_t matrix_size{ rows * columns };
        for (size_t t{ 0 }, k{ 0 }; t < matrix_size; t += rows + 1) {
            temp.values[t] = 1;
        }
        return temp;
    }*/

    friend std::ostream& operator<<(std::ostream& stream, const Matrix& matrix) {
        for (size_t i{ 0 }; i < matrix.values.size(); ++i) {
            stream << matrix.values[i] << ' ';
            if ((i + 1) % matrix.columns == 0)
                stream << '\n';
        }
        return stream;
    }
};


int main()
{
    //SLR
    /*std::vector<Point2D> points{ ReadFile("index.txt", '\t', 2, 3) };
    FindCoeficient(points);
    std::cout << FindCoeficient(points) << std::endl;*/
    try {
        Matrix leftMatrix{ { 1, 2, 3 }, { 4, 5, 6 } };
        Matrix rightMatrix{ { 1, 2}, {3, 4,} };
        std::cout << rightMatrix.Transpose();
        Matrix testmatrix{ {1, 4}, {2, 6} };
        //std::cout << rightMatrix.Determinant();
        //std::cout << testmatrix.Determinant();
       // std::cout << leftMatrix * rightMatrix;
    }
    catch (std::exception& error) {
        std::cout << error.what();
    }
    
    
    Matrix rightMatrix{ { 1, 2, 3}, {4, 12, 6}, {7, 1, 9} };
   // std::cout << rightMatrix.Invertion();
    //MLR
    /*Matrix X("iq_index.txt", '\t', { 2, 3, 4 });
    Matrix Y("iq_index.txt", '\t', { 1 });
    Matrix X_Transposed{ X.Transpose() * (1\n) };
    X = X_Transposed * X;
    X.Invert();
    X = X_Transposed * Y;*/

    // Point2DUnitTest();
}


//container for matrix
//methods for add to matrix
//multiplcaiton matrix
//represent vector as a metrix with 48 rows and 3 column
// matrix with coeficients with 3 rows and 1 column

// linear equations:
//invert of the matrix
//matrix minors
//matrix maltiplication
//determinant



void Point2DUnitTest() {

    double num{ 2.6 };
    double num1{ 3.1 };

    Point2D p1{ num, num1 };
    std::cout << p1.X() << " " << p1.Y() << std::endl;

    Point2D p2(num, num1);
    std::cout << p2.X() << " " << p2.Y() << std::endl;

    Point2D p3 = p1;
    std::cout << p3.X() << " " << p3.Y() << std::endl;

    try {
        Point2D p4{ num, num1, num };
    }
    catch (std::invalid_argument& error) {
        std::cout << error.what();
    }
}

std::vector<Point2D> ReadFile(std::string filename, char delimeter, size_t X_col_index, size_t Y_col_index) {

    if (X_col_index == Y_col_index)
        throw std::invalid_argument("Write X and Y which are not equal.");


    std::ifstream file(filename);
    if (!file.is_open())
        throw std::invalid_argument(std::string("File ") + filename + std::string(" can't open."));

    std::vector<Point2D> Points;

    std::string str;
    std::string currentSmaller;
    std::string currentGreater;
    size_t smallerPosition;
    size_t greaterPosition;


    // from file to string 
    std::getline(file, str, '\n');
    while (!file.eof()) {
        std::getline(file, str, '\n');
        if (str.size() == 0)
            continue;

        std::stringstream ss;
        ss << str;

        size_t NumberOfColumns{ 0 };
        for (size_t i{ 0 }; i < str.size(); ++i) {
            if (str[i] == '\t')
                ++NumberOfColumns;
        }
        ++NumberOfColumns;

        if (X_col_index > NumberOfColumns || Y_col_index > NumberOfColumns) {
            throw std::invalid_argument("Write smaller value");
        }

        if (X_col_index < Y_col_index) {
            smallerPosition = X_col_index;
            greaterPosition = Y_col_index;

            size_t i{ 0 };
            for (; i <= smallerPosition; ++i, std::getline(ss, currentSmaller, delimeter));
            for (size_t j{ i }; j <= greaterPosition; ++j, std::getline(ss, currentGreater, delimeter));

            Points.push_back(Point2D(atof(currentSmaller.c_str()), atof(currentGreater.c_str())));
        }
        else {
            smallerPosition = Y_col_index;
            greaterPosition = X_col_index;

            size_t i{ 0 };
            for (; i <= smallerPosition; ++i, std::getline(ss, currentSmaller, delimeter));
            for (size_t j{ i }; j <= greaterPosition; ++j, std::getline(ss, currentGreater, delimeter));

            Points.push_back(Point2D(atof(currentGreater.c_str()), atof(currentSmaller.c_str())));
        }
    }
    return Points;
}


std::pair<double, double> FindCoeficient(std::vector<Point2D> points) {
    double beta_1;
    double beta_0;

    double X_average{ 0 };
    for (size_t i{ 0 }; i < points.size(); ++i) {
        X_average += points[i].X();
    }
    X_average /= points.size();


    double Y_average{ 0 };
    for (size_t i{ 0 }; i < points.size(); ++i) {
        Y_average += points[i].Y();
    }
    Y_average /= points.size();


    auto CurNumeratorPart = [&X_average, &Y_average](const double& sum, const Point2D& cur)
    {

        return sum + (cur.X() - X_average) * (cur.Y() - Y_average);
    };
    double numerator{ std::accumulate(points.begin(), points.end(), 0.0, CurNumeratorPart) };

    double denominator{ std::accumulate(points.begin(), points.end(), 0.0,
        [&X_average, &Y_average](const double& sum, const Point2D& cur) { return sum + (pow(cur.X() - X_average, 2)); }) };


    beta_1 = numerator / denominator; // 1.373
    beta_0 = Y_average - beta_1 * X_average; // 4.267

    return std::pair<double, double>(beta_1, beta_0);
}





