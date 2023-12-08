//  matrix2d.h
//  Template class matrix2d, two dimensional matrix

#pragma once

template <class T> class matrix2d {
private:
	int nn;   //  number of rows
	int mm;   //  number of columns
	T **v;    //  double pointer to array
public:
	matrix2d();
	matrix2d(int n, int m);			      // Zero-based array
	matrix2d(int n, int m, const T& a);	// Initialize array to constant a
	matrix2d(int n, int m, const T* a);	// Initialize to array a
	matrix2d(const matrix2d &rhs);		// Copy constructor
	matrix2d & operator=(const matrix2d& rhs);	// assignment
	matrix2d & operator=(const T& a);		// assign a to every element
	inline T* operator[](const int i);	   // subscripting: pointer to row i
	inline const T* operator[](const int i) const; // subscripting: pointer to row i
	inline int nrows() const;           //  get number of rows
	inline int ncols() const;           //  get number of columns
   void assign(int,int,const T& a);    //  initialize to constant a
   T sum(int i);                       //  calculate sum over row i
   void SetRows(const int i, const T& a);  //  set row i to a constant a
   void SetCols(const int i, const T& a);  //  set column i to a constant a
	~matrix2d();                        //  destructor
};

template <class T>   //  empty constructor
matrix2d<T>::matrix2d() : nn(0), mm(0), v(0) {}

template <class T>
matrix2d<T>::matrix2d(int n, int m) : nn(n), mm(m), v(new T*[n])
{
	v[0] = new T[m*n];
	for (int i=1; i< n; ++i)
		v[i] = v[i-1] + m;
}

template <class T>
matrix2d<T>::matrix2d(int n, int m, const T& a) : nn(n), mm(m), v(new T*[n])
{
	int i,j;
	v[0] = new T[m*n];
	for (i=1; i< n; ++i)
		v[i] = v[i-1] + m;
	for (i=0; i< n; ++i)
		for (j=0; j<m; ++j)
			v[i][j] = a;
}

template <class T>
matrix2d<T>::matrix2d(int n, int m, const T* a) : nn(n), mm(m), v(new T*[n])
{
	int i,j;
	v[0] = new T[m*n];
	for (i=1; i< n; ++i)
		v[i] = v[i-1] + m;
	for (i=0; i< n; ++i)
		for (j=0; j<m; ++j)
			v[i][j] = *a++;
}

template <class T>
matrix2d<T>::matrix2d(const matrix2d &rhs) : nn(rhs.nn), mm(rhs.mm), v(new T*[nn])
{
	int i,j;
	v[0] = new T[mm*nn];
	for (i=1; i< nn; ++i)
		v[i] = v[i-1] + mm;
	for (i=0; i< nn; ++i)
		for (j=0; j<mm; ++j)
			v[i][j] = rhs[i][j];
}

template <class T>
void matrix2d<T>::assign(int n, int m, const T& a)
{
   nn = n;
   mm = m;
   v = new T*[n];
   register int i,j;
   v[0] = new T[m*n];
	for (i=1; i< n; ++i)
		v[i] = v[i-1] + m;
	for (i=0; i< n; ++i)
		for (j=0; j<m; ++j)
			v[i][j] = a;
}

template <class T>
matrix2d<T> & matrix2d<T>::operator=(const matrix2d<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
	if (this != &rhs) {
		int i,j;
		if (nn != rhs.nn || mm != rhs.mm) {
			if (v != 0) {
				delete[] (v[0]);
				delete[] (v);
			}
			nn=rhs.nn;
			mm=rhs.mm;
			v = new T*[nn];
			v[0] = new T[mm*nn];
		}
		for (i=1; i< nn; ++i)
			v[i] = v[i-1] + mm;
		for (i=0; i< nn; ++i)
			for (j=0; j<mm; ++j)
				v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
matrix2d<T> & matrix2d<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i< nn; ++i)
		for (int j=0; j<mm; ++j)
			v[i][j] = a;
	return *this;
}

template <class T>
inline T* matrix2d<T>::operator[](const int i)	//subscripting: pointer to row i
{
	return v[i];
}

template <class T>
inline const T* matrix2d<T>::operator[](const int i) const // subscripting: pointer to row i
{
	return v[i];
}

template <class T>
T matrix2d<T>::sum(int i)
{
	T* x = v[i];
   T sum = 0;
   for (register int j=0; j<mm; ++j)
      sum += x[j];
   return sum;
}

template <class T>
inline int matrix2d<T>::nrows() const
{
	return nn;
}

template <class T>
inline int matrix2d<T>::ncols() const
{
	return mm;
}

template <class T>
void matrix2d<T>::SetRows(int i, const T& a)
{
   for (register int j=0; j<mm; ++j)
      v[i][j] = a;
}

template <class T>
void matrix2d<T>::SetCols(int i, const T& a)
{
   for (register int j=0; j<mm; ++j)
      v[j][i] = a;
}

template <class T>
matrix2d<T>::~matrix2d()
{
	if (v != 0) {
		delete[] (v[0]);
		delete[] (v);
	}
}
