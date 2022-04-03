#include <delaunator.hpp>
#include <cstdio>
#include <vector>
#include <Eigen/Core>
#include <chrono>
#include "matplotlibcpp.h"
using namespace std;
namespace plt = matplotlibcpp;
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

void add_edge(vector<pair<int, int>>& edges, int p1_index, int p2_index)
{
    vector<pair<int, int>>::const_iterator need_rease;
    bool need_rease_flag = false;
    vector<pair<int, int>>::const_iterator iter;
    for (iter = edges.begin(); iter != edges.end(); iter++)
    {
        bool f1 = iter->first == p1_index && iter->second == p2_index;
        bool f2 = iter->first == p2_index && iter->second == p1_index;
        if (f1 || f2)
        {
            need_rease_flag = true;
            need_rease = iter;
            break;
        }
    }
    if(need_rease_flag)
        edges.erase(need_rease);
    else
        edges.emplace_back(make_pair(p1_index, p2_index));    
}

void add_edge(vector<pair<Eigen::Vector2d, Eigen::Vector2d>>& edges, Eigen::Vector2d & p1, Eigen::Vector2d & p2)
{
    vector<pair<Eigen::Vector2d, Eigen::Vector2d>>::const_iterator need_rease;
    bool need_rease_flag = false;
    vector<pair<Eigen::Vector2d, Eigen::Vector2d>>::const_iterator iter;
    for (iter = edges.begin(); iter != edges.end(); iter++)
    {
        bool f1 = iter->first == p1 && iter->second == p2;
        bool f2 = iter->first == p2 && iter->second == p1;
        if (f1 || f2)
        {
            need_rease_flag = true;
            need_rease = iter;
            break;
        }
    }
    if(need_rease_flag)
        edges.erase(need_rease);
    else
        edges.emplace_back(make_pair(p1, p2));    
}

int main() {
    /* x0, y0, x1, y1, ... */
    string file ="../in.txt";
	//txt文件中有4列
	int columns = 3;
	vector<Eigen::Vector3d> output_vector;
	if (!read_scanf(file, columns, output_vector))
	{
		return 0;
	}
	std::vector<double> coords;
    const auto start = std::chrono::high_resolution_clock::now();
	
    coords.reserve(output_vector.size() * 2);
    // vector<double> x, y;
    // x.reserve(output_vector.size());
    // y.reserve(output_vector.size());
	for (auto & iter : output_vector)
	{
		// points.emplace_back(dt::Vector2<double>{iter(0), iter(1)});
        coords.emplace_back(iter(0));
        coords.emplace_back(iter(1));
        // x.emplace_back(iter(0));
        // y.emplace_back(iter(1));
	}
    // for (auto & iter : coor)
    // {
    //     /* code */
    // }
    // plt::scatter(x, y);
    // plt::show();
    // std::vector<double> coords = {-1, 1, 1, 1, 1, -1, -1, -1};
    double alpha = 0.08;
    //triangulation happens here
    delaunator::Delaunator d(coords);
    cout<<"coords size is "<<coords.size()<<endl;
    auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff = end - start;
    cout<<"cost Delaunator "<<diff.count()<<endl;
    // vector<pair<Eigen::Vector2d, Eigen::Vector2d>> edges;
    vector<pair<int, int>> edges;
    // d.triangles[i] 对应
    for (size_t i = 0; i < d.triangles.size(); i += 3)
    {   
        int p1_index = d.triangles[i];
        int p2_index = d.triangles[i + 1];
        int p3_index = d.triangles[i + 2];

        Eigen::Vector2d p1(d.coords[2 * p1_index], d.coords[2 * p1_index + 1]);
        Eigen::Vector2d p2(d.coords[2 * p2_index], d.coords[2 * p2_index + 1]);
        Eigen::Vector2d p3(d.coords[2 * p3_index], d.coords[2 * p3_index + 1]);
        double a = (p1 - p2).norm();
        double b = (p2 - p3).norm();
        double c = (p3 - p1).norm();
        double s = (a + b + c)/2.0;
        double area = sqrt(s * (s - a) * (s - b) * (s - c));
        double circum_r = a*b*c/(4.0*area);
        if (circum_r < alpha)
        {
            add_edge(edges, p1_index, p2_index);
            add_edge(edges, p2_index, p3_index);
            add_edge(edges, p3_index, p1_index);
            // add_edge(edges, p1, p2);
            // add_edge(edges, p2, p3);
            // add_edge(edges, p3, p1);
        }
    }
    std::cout << edges.size() << " triangles generated in " << diff.count()
			<< "s\n";
    end = std::chrono::high_resolution_clock::now();
	diff = end - start;
    // vector<double> z1,w1;
    // z1.reserve(2);
    // w1.reserve(2);
    // cout<<"from "<<edges.front().first.transpose()<<" to "<<edges.front().second.transpose()<<endl;
    // z1.emplace_back(edges.front().first(0));
    // z1.emplace_back(edges.front().second(0));
    // w1.emplace_back(edges.front().first(1));
    // w1.emplace_back(edges.front().second(1));
    // plt::plot(z1,w1);

    // vector<double> z2,w2;
    // z2.reserve(2);
    // w2.reserve(2);
    // z2.emplace_back(edges.at(1).first(0));
    // z2.emplace_back(edges.at(1).second(0));
    // w2.emplace_back(edges.at(1).first(1));
    // w2.emplace_back(edges.at(1).second(1));
    // plt::plot(z2,w2);  
    // vector<double> z,w;
    // for (auto & iter : edges)
    // {
    //     z.clear();
    //     w.clear();
    //     z.emplace_back(iter.first(0));
    //     z.emplace_back(iter.second(0));
    //     w.emplace_back(iter.first(1));
    //     w.emplace_back(iter.second(1));
    //     plt::plot(z,w);
    // }
    
    // plt::show();
	
    cout<<"outliner size is "<<edges.size()<<endl;
    // for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
    //     printf(
    //         "Triangle points: [[%f, %f], [%f, %f], [%f, %f]]\n",
    //         d.coords[2 * d.triangles[i]],        //tx0
    //         d.coords[2 * d.triangles[i] + 1],    //ty0
    //         d.coords[2 * d.triangles[i + 1]],    //tx1
    //         d.coords[2 * d.triangles[i + 1] + 1],//ty1
    //         d.coords[2 * d.triangles[i + 2]],    //tx2
    //         d.coords[2 * d.triangles[i + 2] + 1] //ty2
    //     );
    // }

}
