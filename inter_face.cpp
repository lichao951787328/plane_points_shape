#include "PythonParam.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <Eigen/Core>
#include <list>
#include <SFML/Graphics.hpp>
#include "vector2.h"
#include "triangle.h"
#include "delaunay.h"
using namespace std;
#define Point Eigen::Vector2d
int read_scanf(const string &filename, const int &cols, vector<Eigen::Vector3d> &_vector)
{
	FILE *fp = fopen(filename.c_str(), "r");
	bool flag = true;
	int i = 0;
	if (!fp) 
	{ 
		cout << "File open error!\n"; 
		return 0; 
	}
	while (flag)
	{
		Eigen::Vector3d rowArray;
		// double *rowArray = new double[cols]; //new一个double类型的动态数组
		for (i = 0; i < cols; i++) //读取数据，存在_vector[cols]中
		{
			if (EOF == fscanf(fp,"%lf", &rowArray(i)))
			{ 
				flag = false; 
				break; 
			}
		}
		if (cols == i) //将txt文本文件中的一行数据存入rowArray中，并将rowArray存入vector中
			_vector.push_back(rowArray);
	}
	fclose(fp);
	return 1;
}

int main(){
	// 点云生成
	string file ="../in.txt";
	//txt文件中有4列
	int columns = 3;
	vector<Eigen::Vector3d> output_vector;
	if (!read_scanf(file, columns, output_vector))
	{
		return 0;
	}
	
	cout<<"number of points "<<output_vector.size()<<endl;
	cout<<output_vector.front().transpose()<<endl;
	cout<<output_vector.back().transpose()<<endl;

	vector<Eigen::Vector2d> points;
	points.reserve(output_vector.size());
	for (auto & iter : output_vector)
	{
		points.emplace_back(iter.head(2));
	}
	cout<<points.front().transpose()<<endl;
	cout<<points.back().transpose()<<endl;

	PythonParamBuilder pointPara;
	pointPara.AddVectorPoint(points);
	PyObject* point_py = pointPara.Build();
	PyObject* alpha_py = Py_BuildValue("d", 0.03);

	PyObject* pArgs = PyTuple_New(2);
    PyTuple_SetItem(pArgs, 0, point_py); 
    PyTuple_SetItem(pArgs, 1, alpha_py); 
	cout<<"has create param "<<endl;
	PyObject* pRet = PyObject_CallObject(pFunc, pArgs);
	cout<<"back to cpp"<<endl;
	int res = 0;
	PyArg_Parse(pRet, "i", &res );
	cout<<"res: "<<res<<endl;
	Py_ssize_t size_ = PyTuple_Size(pRet);
	cout<<"size_ is "<<size_<<endl;
    // Py_DECREF(pArgs);
	Py_DECREF(pModule);
	Py_DECREF(pFunc);
    Py_DECREF(pRet);
	cout<<"back to cpp"<<endl;

	// int res = 0;
	// PyArg_Parse(pRet, "i", &res );
    //     Py_DECREF(pRet);
	// cout << res << endl;

	// Py_Finalize();

	// PyObject_CallObject(pFunc, NULL );

	Py_Finalize();
	cout<<"back to cpp"<<endl;

	return 0;
}
