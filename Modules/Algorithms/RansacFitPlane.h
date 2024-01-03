#pragma once
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <functional>
#include <random>
#include <Eigen/Eigen>
using namespace std;

struct float3
{
    float x, y, z;
};

struct plane
{
    float x, y, z, D;
};

//不建议使用该方法，因为下面有更完善的方法
void ransacFitPlane2(std::vector<Eigen::Vector3d>& pts_3d, int max_iter, float threshold)
{
    srand(time(0)); //随机种子
    int size_old = 3;
    double a, b, c, d; //平面法向量系数
    double outa, outb, outc, outd; //平面法向量系数

    while (--max_iter) //设置循环的次数
    {
        vector<int> index;
        for (int k = 0; k < 3; ++k)
        {
            index.push_back(rand() % pts_3d.size()); //随机选取三个点 
        }
        auto idx = index.begin();
        double x1 = pts_3d.at(*idx)[0], y1 = pts_3d.at(*idx)[1], z1 = pts_3d.at(*idx)[2];
        ++idx;
        double x2 = pts_3d.at(*idx)[0], y2 = pts_3d.at(*idx)[1], z2 = pts_3d.at(*idx)[2];
        ++idx;
        double x3 = pts_3d.at(*idx)[0], y3 = pts_3d.at(*idx)[1], z3 = pts_3d.at(*idx)[2];

        a = (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1);
        b = (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1);
        c = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
        d = -(a * x1 + b * y1 + c * z1);

        for (auto iter = pts_3d.begin(); iter != pts_3d.end(); ++iter)
        {
            double dis = fabs(a * iter[0] + b * iter[1] + c * iter[2] + d) / sqrt(a * a + b * b + c * c);//点到平面的距离公式
            if (dis < threshold) index.push_back(iter - pts_3d.begin());
        }
        //更新集合
        if (index.size() > size_old)
        {
            size_old = index.size();
            outa = a;
            outb = b;
            outc = c;
            outd = d;
        }
        index.clear();
    }
    cout << outa << " " << outb << " " << outc << " " << outd << endl;
}

inline void ImportPointCloud(char& filename, std::vector<float3>& points)
{
    float3 temp;
    bool mStatus = false;
    std::ifstream import;
    import.open(&filename);
    while (import.peek() != EOF)
    {
        import >> temp.x >> temp.y >> temp.z;
        points.push_back(temp);
    }
    import.close();
}

inline float DistanceToPlane(float3& point, plane& mPlane)
{
    return std::abs(mPlane.x * point.x + mPlane.y * point.y + mPlane.z * point.z + mPlane.D);
}

inline bool PlaneFromPoints(std::vector<float3>& points, plane& mPlane)
{
    if (points.size() < 3) {
        throw std::runtime_error("Not enough points to calculate plane");
        return false;
    }

    float3 sum = { 0,0,0 };

    for (int i = 0; i < static_cast<int>(points.size()); i++)
    {
        sum.x += points[i].x;
        sum.y += points[i].y;
        sum.z += points[i].z;
    }

    float3 centroid = { 0,0,0 };
    centroid.x = sum.x / float(points.size());
    centroid.y = sum.y / float(points.size());
    centroid.z = sum.z / float(points.size());


    float xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
    for (int i = 0; i < static_cast<int>(points.size()); i++)
    {
        float3 temp;
        temp.x = points[i].x - centroid.x;
        temp.y = points[i].y - centroid.y;
        temp.z = points[i].z - centroid.z;

        xx += temp.x * temp.x;
        xy += temp.x * temp.y;
        xz += temp.x * temp.z;
        yy += temp.y * temp.y;
        yz += temp.y * temp.z;
        zz += temp.z * temp.z;
    }

    float detX = yy * zz - yz * yz;
    float detY = xx * zz - xz * xz;
    float detZ = xx * yy - xy * xy;

    float detMax = std::max(std::max(detX, detY), detZ);
    if (detMax <= 0)
    {
        mPlane.x = 0; mPlane.y = 0; mPlane.z = 0; mPlane.D = 0;
        return false;
    }

    float3 dir{};
    if (detMax == detX)
    {
        float a = static_cast<float>((xz * yz - xy * zz) / detX);
        float b = static_cast<float>((xy * yz - xz * yy) / detX);
        dir = { 1.0, a, b };
    }
    else if (detMax == detY)
    {
        float a = static_cast<float>((yz * xz - xy * zz) / detY);
        float b = static_cast<float>((xy * xz - yz * xx) / detY);
        dir = { a, 1.0, b };
    }
    else
    {
        float a = static_cast<float>((yz * xy - xz * yy) / detZ);
        float b = static_cast<float>((xz * xy - yz * xx) / detZ);
        dir = { a, b, 1.0 };
    }

    float dirTemp = sqrt(dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
    dir.x /= dirTemp;
    dir.y /= dirTemp;
    dir.z /= dirTemp;

    mPlane.x = dir.x;
    mPlane.y = dir.y;
    mPlane.z = dir.z;
    mPlane.D = -(mPlane.x * centroid.x + mPlane.y * centroid.y + mPlane.z * centroid.z);
    return true;
}

void Ransac3DPlane(std::vector<float3>& inPoints, plane& outPlane)
{
    float minSamplePercent = 0.4f;  // Select the initial sample percentage
    int minSampleNum = static_cast<int>(float(inPoints.size()) * minSamplePercent);
    float maxSamplePercent = 0.9f;  // Maximum percentage control iteration stop
    int maxSampleNum = static_cast<int>(float(inPoints.size()) * maxSamplePercent);

    float distThreshold = 0.2f;
    int iterations = 0;
    int maxIterations = 30;
    std::vector<float3> maxPlanePoints; // Satisfy the maximum number of interior points of the plane model

    while (iterations < maxIterations || maxPlanePoints.size() < maxSampleNum)
    {
        iterations++;
        float3 ptTemp1, ptTemp2, ptTemp3;
        plane planeParame = { 0, 0, 0, 0 };
        std::vector<float3> points;
        std::vector<float3> tempPlanePoints;
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(inPoints.begin(), inPoints.end(),g);
        for (int i = 0; i < minSampleNum; i++)
        {
            ptTemp1.x = inPoints[i].x;
            ptTemp1.y = inPoints[i].y;
            ptTemp1.z = inPoints[i].z;

            points.push_back(ptTemp1);
        }

        if(!PlaneFromPoints(points, planeParame))
            continue;
        points.clear();

        for (int i = 0; i < static_cast<int>(inPoints.size()); i++)
        {
            float distTemp;
            ptTemp2.x = inPoints[i].x;
            ptTemp2.y = inPoints[i].y;
            ptTemp2.z = inPoints[i].z;
            distTemp = DistanceToPlane(ptTemp2, planeParame);
            if (distTemp < distThreshold)
            {
                tempPlanePoints.push_back(ptTemp2);
            }
            else { continue; }
        }

        if (tempPlanePoints.size() > maxPlanePoints.size())
        {
            maxPlanePoints.clear();
            maxPlanePoints.reserve(tempPlanePoints.size());
            for (int i = 0; i < static_cast<int>(tempPlanePoints.size()); i++)
            {
                ptTemp3.x = tempPlanePoints[i].x;
                ptTemp3.y = tempPlanePoints[i].y;
                ptTemp3.z = tempPlanePoints[i].z;
                maxPlanePoints.push_back(ptTemp3);
            }
        }
        tempPlanePoints.clear();

    }

    PlaneFromPoints(maxPlanePoints, outPlane);
}

void test()
{
    char filename[255];
    sprintf_s(filename, 255, "111.asc");
    std::vector<float3> points;
    plane mPlane;
    ImportPointCloud(*filename, points);
    //PlaneFromPoints(points, mPlane);
    Ransac3DPlane(points, mPlane);
    std::cout << mPlane.x << "  " << mPlane.y << "  " << mPlane.z << "  " << mPlane.D << std::endl;

}

//int main()
//{
//    test();
//    return 0;
//}
