/*
Copied from 

https://github.com/OmarAflak/Matrix/blob/master/matrix.h

---------------------------------------------------------

Usage:

#include "matrix.h"
...

// Constructors :
Matrix<type_of_variable> matrix(int height, int width);
Matrix<type_of_variable> matrix(vector<vector<type_of_variable> > &array);
Matrix<type_of_variable> matrix();

Operators

Matrix<int> a,b,c;
int z;

c = a + b;
c = a - b;
c = a * b;
c = z * a;
c = a / z;

Example

Matrix<int> A(4,5);
Matrix<int> B(4,5);
Matrix<int> C(5,6);

Matrix<int> D = A + B;  // = A.add(B)
Matrix<int> D = A - B;  // = A.subtract(B)
Matrix<int> D = A * B;  // = A.multiply(B)
Matrix<int> D = B.dot(C);
Matrix<int> D = A.transpose();
Matrix<int> D = A.subMatrix(0,0, 2,2); // get a sub matrix of size 2x2 from coordinates 0;0 in A

A += B;
A -= B;
A *= B;

int multiplyByTwo(int x){
    return 2*x;
}

Matrix<int> D = A.applyFunction(multiplyByTwo);

// you can combine operations :

Matrix<int> D = (A+B).dot(C).applyFunction(multiplyByTwo);

std::cout << D << std::endl;
----------------------------------------------------------
*/

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>

template <class data_t>
class Matrix{
    private:
        std::vector<std::vector<data_t>> array; // 2D vector type to store the matrix values
        int height = 2;
        int width = 2;

    public:
	// Make nan
	void MakeNaN(){
		for(int j=0; j<height; j++){
			for(int k=0; k<width; k++){
				array[j][k] = nan;
			}
		}
	}

	// find the condition number for inversion for 2x2 matrices
	double cond(){
		data_t S1 = pow(array[0][0], 2) + pow(array[1][0], 2) + pow(array[0][1], 2) + pow(array[1][1], 2);
		data_t acbd = array[0][0]*array[1][0] + array[0][1]*array[1][1];
		data_t S2 = std::sqrt(pow(pow(array[0][0], 2) + pow(array[0][1], 2) - pow(array[1][0], 2) - pow(array[1][1], 2),2) + 4*pow(acbd, 2));
		double sigma1 = std::sqrt(std::abs(S1 + S2));
		double sigma2 = std::sqrt(std::abs(S1 - S2));
		double sigma_min = std::min(sigma1, sigma2);
		if (sigma_min == 0)
		{
			return INFINITY;
		} else {
			double output = std::max(sigma1, sigma2)/sigma_min;
			return output;
		}
	}

	// sum over all indices
	data_t sum(){
		data_t sum_ = 0;
		for(int j=0; j<height; j++){
                        for(int k=0; k<width; k++){
                                sum_ += array[j][k];
                        }
                }
		return sum_;
	}

	// compute trace
        data_t trace(){
		if (width==height){
                	data_t trace_ = 0;
	                for(int j=0; j<height; j++){
	                        sum_ += array[j][j];
	                }
	                return trace_;
		} else {
			throw std::invalid_argument("trace() only works on square matrices");
        	}
	}
	
	// Initialise 2x2 matrix
	void init(data_t a00, data_t a01, data_t a10, data_t a11){
		array[0][0] = a00;
		array[1][0] = a10;
		array[0][1] = a01;
		array[1][1] = a11;
	}

	/*
	// Initialise general matrix from a width x height array
	void init(data_t init_arr[width][height]){
		for(int j=0; j<height; j++){
                        for(int k=0; k<width; k++){
                                array[j][k] = init_arr[j][k];
                        }
                }
        }*/		

	// 2x2 Inverse
	Matrix<data_t> inverse(){
		if ( height==2 && width ==2 ){
			Matrix<data_t> Inv(2, 2);
			data_t determinant = array[0][0]*array[1][1] - array[0][1]*array[1][0];
			if (abs(determinant) == 0){
				throw std::overflow_error("Inverse(): determinant = zero exception");
			}
			else {
				Inv.array[0][0] = array[1][1]/determinant;
				Inv.array[1][0] = -array[1][0]/determinant;
				Inv.array[0][1] = -array[0][1]/determinant;
				Inv.array[1][1] = array[0][0]/determinant;
				return Inv;
			}
		}
		else {
			throw std::invalid_argument("Inverse() only works on 2x2 matrices");
		}
	}

	// -- copied code -- //

	// three different constructor options
        Matrix<data_t>(int height, int width);
        Matrix<data_t>(std::vector<std::vector<data_t> > const &array);
        Matrix<data_t>();

        int getHeight() const;
        int getWidth() const;

        Matrix<data_t> add(const Matrix<data_t>& m) const;
        Matrix<data_t> subtract(const Matrix<data_t>& m) const;
        Matrix<data_t> multiply(const Matrix<data_t>& m) const;
        Matrix<data_t> dot(const Matrix<data_t>& m) const;
        Matrix<data_t> transpose() const;
        Matrix<data_t> multiply(const data_t& value) const;
        Matrix<data_t> divide(const data_t& value) const;

        Matrix<data_t> applyFunction(data_t (*function)(data_t)) const;
        Matrix<data_t> subMat(int startH, int startW, int h, int w) const;

        void fill(const data_t& value);
        void put(int h, int w, const data_t& value);
        data_t get(int h, int w) const;

        void print(std::ostream &flux) const;

        bool operator==(const Matrix<data_t>& m);
        bool operator!=(const Matrix<data_t>& m);
        Matrix<data_t> operator+=(const Matrix<data_t>& m);
        Matrix<data_t> operator-=(const Matrix<data_t>& m);
        Matrix<data_t> operator*=(const Matrix<data_t>& m);
        Matrix<data_t> operator*=(const data_t &m);
        Matrix<data_t> operator/=(const data_t &m);
        data_t& operator()(int y, int x);
};

template <class data_t> Matrix<data_t> operator+(const Matrix<data_t>& a, const Matrix<data_t>& b);
template <class data_t> Matrix<data_t> operator-(const Matrix<data_t>& a, const Matrix<data_t>& b);
template <class data_t> Matrix<data_t> operator*(const Matrix<data_t>& a, const Matrix<data_t>& b);
template <class data_t> Matrix<data_t> operator*(const data_t &b, const Matrix<data_t>& a);
template <class data_t> Matrix<data_t> operator/(const Matrix<data_t>& a, const data_t &b);
template <class data_t> std::ostream& operator<<(std::ostream &flux, const Matrix<data_t>& m);

template <class data_t>
Matrix<data_t>::Matrix(int height, int width){
    this->height = height;
    this->width = width;
    this->array = std::vector<std::vector<data_t> >(height, std::vector<data_t>(width));
}

template <class data_t>
Matrix<data_t>::Matrix(std::vector<std::vector<data_t> > const &array){
    if(array.size()==0)
        throw std::invalid_argument("Size of array must be greater than 0.");

    this->height = array.size();
    this->width = array[0].size();
    this->array = array;
}

template <class data_t>
Matrix<data_t>::Matrix(){
    height = 0;
    width = 0;
}

template <class data_t>
int Matrix<data_t>::getHeight() const{
    return height;
}

template <class data_t>
int Matrix<data_t>::getWidth() const{
    return width;
}

template <class data_t>
void Matrix<data_t>::fill(const data_t& value){
    for (int i=0 ; i<height ; i++){
        for (int j=0 ; j<width ; j++){
            array[i][j] = value;
        }
    }
}

template <class data_t>
void Matrix<data_t>::put(int h, int w, const data_t& value){
    if(!(h>=0 && h<height && w>=0 && w<width))
        throw std::invalid_argument("Index out of bounds.");

    array[h][w] = value;
}

template <class data_t>
data_t Matrix<data_t>::get(int h, int w) const{
    if(!(h>=0 && h<height && w>=0 && w<width))
        throw std::invalid_argument("Index out of bounds.");

    return array[h][w];
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::multiply(const data_t& value) const{
    Matrix result(array);
    for (int i=0 ; i<height ; i++){
        for (int j=0 ; j<width ; j++){
            result.array[i][j] *= value;
        }
    }

    return result;
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::divide(const data_t& value) const{
    Matrix result(array);
    for (int i=0 ; i<height ; i++){
        for (int j=0 ; j<width ; j++){
            result.array[i][j] /= value;
        }
    }

    return result;
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::add(const Matrix& m) const{
    if(!(height==m.height && width==m.width))
        throw std::invalid_argument("Matrix dimension must be the same.");

    Matrix result(height, width);
    for (int i=0 ; i<height ; i++){
        for (int j=0 ; j<width ; j++){
            result.array[i][j] = array[i][j] + m.array[i][j];
        }
    }

    return result;
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::subtract(const Matrix& m) const{
    if(!(height==m.height && width==m.width))
        throw std::invalid_argument("Matrix dimension must be the same.");

    Matrix result(height, width);
    for (int i=0 ; i<height ; i++){
        for (int j=0 ; j<width ; j++){
            result.array[i][j] = array[i][j] - m.array[i][j];
        }
    }
    return result;
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::multiply(const Matrix& m) const{
    if(!(height==m.height && width==m.width))
        throw std::invalid_argument("Matrix dimension must be the same.");

    Matrix result(height, width);

    for (int i=0 ; i<height ; i++){
        for (int j=0 ; j<width ; j++){
            result.array[i][j] = array[i][j] * m.array[i][j];
        }
    }
    return result;
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::dot(const Matrix& m) const{
    if(!(width==m.height))
        throw std::invalid_argument("Dot product not compatible.");

    data_t w=0;
    int mwidth = m.width;

    Matrix<data_t> result(height, mwidth);
    for (int i=0 ; i<height ; i++){
        for (int j=0 ; j<mwidth ; j++){
            for (int h=0 ; h<width ; h++){
                w += array[i][h]*m.array[h][j];
            }
            result.array[i][j] = w;
            w=0;
        }
    }

    return result;
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::transpose() const{
    Matrix<data_t> result(width, height);

    for (int i=0 ; i<width ; i++){
        for (int j=0 ; j<height ; j++){
            result.array[i][j] = array[j][i];
        }
    }
    return result;
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::applyFunction(data_t (*function)(data_t)) const{
    Matrix<data_t> result(height, width);

    for (int i=0 ; i<height ; i++){
        for (int j=0 ; j<width ; j++){
            result.array[i][j] = (*function)(array[i][j]);
        }
    }

    return result;
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::subMat(int startH, int startW, int h, int w) const{
    if(!(startH>=0 && startH+h<=height && startW>=0 && startW+w<=width))
        throw std::invalid_argument("Index out of bounds");

    Matrix<data_t> result(h,w);
    for (int i=startH ; i<startH+h ; i++){
        for (int j=startW ; j<startW+w ; j++){
            result.array[i-startH][j-startW] = array[i][j];
        }
    }
    return result;
}

template <class data_t>
void Matrix<data_t>::print(std::ostream &flux) const{
    int maxLength[width] = {};
    std::stringstream ss;

    for (int i=0 ; i<height ; i++){
        for (int j=0 ; j<width ; j++){
            ss << array[i][j];
            if(maxLength[j] < ss.str().size()){
                maxLength[j] = ss.str().size();
            }
            ss.str(std::string());
        }
    }

    for (int i=0 ; i<height ; i++){
        for (int j=0 ; j<width ; j++){
            flux << array[i][j];
            ss << array[i][j];
            for (int k=0 ; k<maxLength[j]-ss.str().size()+1 ; k++){
                flux << " ";
            }
            ss.str(std::string());
        }
        flux << std::endl;
    }
}

template <class data_t>
bool Matrix<data_t>::operator==(const Matrix& m){
    if(height==m.height && width==m.width){
        for (int i=0 ; i<height ; i++){
            for (int j=0 ; j<width ; j++){
                if(array[i][j]!=m.array[i][j]){
                    return false;
                }
            }
        }
        return true;
    }
    return false;
}

template <class data_t>
bool Matrix<data_t>::operator!=(const Matrix& m){
    return !operator==(m);
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::operator+=(const Matrix& m){
    this->array = add(m).array;
    return *this;
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::operator-=(const Matrix& m){
    this->array = subtract(m).array;
    return *this;
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::operator*=(const Matrix& m){
    this->array = multiply(m).array;
    return *this;
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::operator*=(const data_t &m){
    *this = this->multiply(m);
    return *this;
}

template <class data_t>
Matrix<data_t> Matrix<data_t>::operator/=(const data_t &m){
    *this = this->divide(m);
    return *this;
}

template <class data_t>
data_t& Matrix<data_t>::operator()(int y, int x){
    if(!(y>=0 && y<height && x>=0 && x<width))
        throw std::invalid_argument("Index out of bounds.");
    return array[y][x];
}

template <class data_t>
Matrix<data_t> operator+(const Matrix<data_t>& a, const Matrix<data_t>& b){
    return a.add(b);
}

template <class data_t>
Matrix<data_t> operator-(const Matrix<data_t>& a, const Matrix<data_t>& b){
    return a.subtract(b);
}

template <class data_t>
Matrix<data_t> operator*(const Matrix<data_t>& a, const Matrix<data_t>& b){
    return a.multiply(b);
}

template <class data_t>
Matrix<data_t> operator*(const data_t &b, const Matrix<data_t>& a){
    return a.multiply(b);
}
template <class data_t>
Matrix<data_t> operator/(const Matrix<data_t>& a, const data_t &b){
    return a.divide(b);
}

template <class data_t>
std::ostream& operator<<(std::ostream &flux, const Matrix<data_t>& m){
    m.print(flux);
    return flux;
}

#endif /* MATRIX_HPP_ */
