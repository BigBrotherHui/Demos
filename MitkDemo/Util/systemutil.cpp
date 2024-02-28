#pragma execution_character_set("utf-8")
#include "systemutil.h"

#include <WinSock2.h>
#include <stdio.h>
#include <ws2tcpip.h>
#include <QApplication>
#include <QDebug>
#include <QDir>
#include <QFile>
#include <QMetaEnum>
#include <QObject>
#include <QTranslator>
#include <cmath>
#include <iostream>
#include <QTime>

//#include "configutil.h"
//#include "logger.h"

#ifndef tr
#define tr QObject::tr
#endif

namespace SystemUtil {
/**
 * @brief get list of child object
 * @param parent
 * @return child list
 */
QObjectList getChild(QObject *parent)
{
    QObjectList list;
    for (QObject *object : parent->children())
    {
        list.append(object);
        if (object->metaObject()->className() == QStringLiteral("QWidget"))
        {
            list.append(getChild(object));
        }
    }
    return list;
}

void InstallEventFitler(QObject *parent, QObjectList objectList)
{
    if (objectList.size() <= 0)
    {
        objectList = getChild(parent);
    }
    for (QObject *object : objectList)
    {
        object->installEventFilter(parent);
    }
}

void LoadStyle() {
//    LOGGER_INFO("Load style");
//  QStringList styleFiles = ConfigUtil::GetStyleFile();
    QString qss;
    QFile qss_file(":/qss/style.css");
    if (qss_file.open(QFile::ReadOnly))
    {
        qss = qss_file.readAll();
    }
//  for (QString styleFile : styleFiles) {
//    QFile qssFile(styleFile);
//    qssFile.open(QFile::ReadOnly);
//    if (qssFile.isOpen()) {
//      LOGGER_INFO("style file loaded:" + styleFile);
//      qss += QLatin1String(qssFile.readAll());
//      qssFile.close();
//    } else {
//      LOGGER_ERROR("style file load failed:" + styleFile);
//    }
//  }
    qApp->setStyleSheet(qss);
}

//void LoadTranslator()
//{
////    LOGGER_INFO("Load translator");
//    QStringList translatorFiles = ConfigUtil::GetTranslatorFile();
//    for (QString translatorFile : translatorFiles)
//    {
//        LOGGER_INFO("Translator file loaded:" + translatorFile);
//        QTranslator *translator = new QTranslator(qApp);
//        if (translator->load(translatorFile))
//        {
//            qApp->installTranslator(translator);
//        }
//        else
//        {
//            LOGGER_ERROR("Translator file load failed:" + translatorFile);
//        }
//    }
//    return;
//}

bool CheckDir(QString pDir, bool pCreate)
{
    return CheckDir(pDir, tr("目录"), pCreate);
}

bool CheckDir(QString pDir, QString pDirDesc, bool pCreate)
{
//    LOGGER_INFO(QString("CheckDir:%1").arg(pDir));
    if (pDirDesc.isEmpty())
    {
        pDirDesc = tr("目录");
    }
    QDir qDir;
    if (qDir.exists(pDir))
    {
        return true;
    }
    if (!pCreate)
    {
        return false;
    }
//    LOGGER_INFO(tr("%1(%2)不存在，自动创建...").arg(pDirDesc).arg(pDir));
    return CreateDir(pDir, pDirDesc);
}

bool CreateDir(QString pDir, QString pDirDesc)
{
    if (pDirDesc.isEmpty())
    {
        pDirDesc = tr("目录");
    }
    QDir qDir;
    if (qDir.exists(pDir))
    {
//        LOGGER_ERROR(tr("%1(%2)已经存在，创建失败").arg(pDirDesc).arg(pDir));
        return false;
    }
    if (!qDir.mkpath(pDir))
    {
//        LOGGER_ERROR(tr("创建%1(%2)失败").arg(pDirDesc).arg(pDir));
        return false;
    }
//    LOGGER_INFO(tr("%1(%2)创建成功").arg(pDirDesc).arg(pDir));
    return true;
}

bool DeleteDir(QString pDir, QString pDirDesc)
{
    qDebug() << "SystemUtil::DeleteDir" << pDir;
    if (pDir.isEmpty())
    {
        return true;
    }
    if (pDirDesc.isEmpty())
    {
        pDirDesc = tr("目录");
    }
    QDir qDir(pDir);
    if (!qDir.exists())
    {
//        LOGGER_INFO(tr("%1(%2)不存在，删除失败").arg(pDirDesc).arg(pDir));
        return true;
    }
    if (!qDir.removeRecursively())
    {
//        LOGGER_ERROR(tr("删除%1(%2)失败").arg(pDirDesc).arg(pDir));
        return false;
    }
//    LOGGER_INFO(tr("%1(%2)删除成功").arg(pDirDesc).arg(pDir));
    return true;
}

QString GetFullDir(QString pDir)
{
    auto dir = pDir.replace("\\", "/");
    if (pDir.indexOf(":") < 0 && !pDir.startsWith("/"))
    {
        //表示为相对路径，则加上当前路径
        auto parentDir = QCoreApplication::applicationDirPath();
        while(dir.startsWith("./"))
        {
            dir = dir.mid(2);
        }
        while(dir.startsWith("../"))
        {
            dir = dir.mid(3);
            parentDir = GetDir(parentDir);
        }
        dir = parentDir + "/" + dir;
    }
    return dir;
}

bool FileCopy(QString pSourceFile, QString pTargetFile, bool pOveride)
{
    if (QFile::exists(pTargetFile))
    {
        if (!pOveride)
        {
//            LOGGER_ERROR(tr("复制文件失败(从%1到%2): 文件已存在")
//                             .arg(pSourceFile)
//                             .arg(pTargetFile));
            return false;
        }
        if (!QFile::remove(pTargetFile))
        {
//            LOGGER_ERROR(tr("复制文件失败(从%1到%2): 删除已有文件失败")
//                             .arg(pSourceFile)
//                             .arg(pTargetFile));
            return false;
        }
    }
    return QFile::copy(pSourceFile, pTargetFile);
}

int CopyDir(QString pSourceDir, QString pTargetDir, bool pOveride)
{
    QStringList filters;
    return CopyDir(pSourceDir, pTargetDir, filters, pOveride);
}

int CopyDir(QString pSourceDir, QString pTargetDir, QString pSuffix, bool pOveride)
{
    QStringList filters;
    filters<<QString("*.%1").arg(pSuffix);
    return CopyDir(pSourceDir, pTargetDir, filters, pOveride);
}

int CopyDir(QString pSourceDir, QString pTargetDir, QStringList pFilters, bool pOveride)
{
    QDir source(pSourceDir);
    if (!source.exists())
    {
//        LOGGER_ERROR(tr("复制目录失败(从%1到%2)").arg(pSourceDir).arg(pTargetDir));
        return -1;
    }
    QDir target(pTargetDir);
    if (target.exists())
    {
        if (!pOveride)
        {
//            LOGGER_ERROR(tr("复制目录失败(从%1到%2)，目标目录已存在")
//                             .arg(pSourceDir)
//                             .arg(pTargetDir));
            return -1;
        }
    }
    else if (!CreateDir(pTargetDir, tr("目标目录")))
    {
        return -1;
    }
    //遍历源文件夹下所有文件
    int fileCnt = 0;
    QStringList files = source.entryList(pFilters, QDir::Files);
    for (QString file : files)
    {
        QString sourceFile = pSourceDir + "/" + file;
        QString targetFile = pTargetDir + "/" + file;
        if(QFile::exists(targetFile))
        {
            QFile::remove(targetFile);
        }
        if (!QFile::copy(sourceFile, targetFile))
        {
//            LOGGER_ERROR(tr("复制文件失败(从%1到%2)")
//                         .arg(sourceFile)
//                         .arg(targetFile));
            return -1;
        }
        fileCnt ++;
    }
    //遍历源文件夹下所有子文件夹
    QStringList filters;
    QStringList dirs = source.entryList(filters, QDir::Dirs | QDir::NoDotAndDotDot);
    for (QString dir : dirs)
    {
        QString subSourceDir = pSourceDir + "/" + dir;
        QString subTargetDir = pTargetDir + "/" + dir;
        int subFileCnt = CopyDir(subSourceDir, subTargetDir, pFilters, pOveride);
        if (subFileCnt==-1) {
          return -1;
        }
        fileCnt += subFileCnt;
    }
    return fileCnt;
}

bool EmptyDir(QString pDir, QString pSuffix, bool pEmptySub)
{
    QStringList filter;
    filter<<QString("*.%1").arg(pSuffix);
    return EmptyDir(pDir, pEmptySub, filter);
}

bool EmptyDir(QString pDir, bool pEmptySub, QStringList pNameFilters)
{
    qDebug() << "EmptyDir" << pDir;
    QDir dir(pDir);
    if (!dir.exists())
    {
        return true;
    }
    //遍历源文件夹下所有文件
    QStringList files = dir.entryList(pNameFilters, QDir::Files);
    for (QString file : files)
    {
        qDebug() << "EmptyDir" << file;
        QFile::remove(QString("%1/%2").arg(pDir).arg(file));
    }
    if (pEmptySub)
    {
        //遍历源文件夹下所有子文件夹
        QStringList dirs = dir.entryList(pNameFilters, QDir::Dirs | QDir::NoDotAndDotDot);
        for (QString dir : dirs)
        {
            EmptyDir(dir, true, pNameFilters);
            if(IsDirEmpty(dir))
            {
                QDir qDir(dir);
                qDir.removeRecursively();
            }
        }
    }
    return true;
}

bool IsDirEmpty(QString pDir)
{
    qDebug() << "IsDirEmpty" << pDir;
    QDir dir(pDir);
    if (!dir.exists())
    {
        return false;
    }
    if (dir.entryList(QDir::Files | QDir::Dirs | QDir::NoDotAndDotDot).size() > 0)
    {
        return false;
    }
    return true;
}

bool CheckLogDir(QString pDir)
{
    QDir qDir;
    if (qDir.exists(pDir))
    {
        return true;
    }
    if (!qDir.mkpath(pDir))
    {
        return false;
    }
    return true;
}

QStringList GetLocalIp()
{
    WORD wVersion;
    WSADATA WSAData;
    wVersion = MAKEWORD(2, 2);
    WSAStartup(wVersion, &WSAData);
    struct addrinfo hints;
    struct addrinfo *res, *cur;
    int ret;
    struct sockaddr_in *addr;
    QStringList ipList;

    memset(&hints, 0, sizeof(struct addrinfo));
    hints.ai_family = AF_INET;   /* Allow IPv4 */
    hints.ai_flags = AI_PASSIVE; /* For wildcard IP address */
    hints.ai_protocol = 0;       /* Any protocol */
    hints.ai_socktype = SOCK_STREAM;
    char hostname[256];
    //获取本机主机名称
    int iRet;
    iRet = gethostname(hostname, 256);
    qDebug() << hostname;
    ret = getaddrinfo(hostname, NULL, &hints, &res);
    qDebug() << ret;

    if (ret != 0)
    {
        perror("getaddrinfo");
        return ipList;
    }
    for (cur = res; cur != NULL; cur = cur->ai_next)
    {
        addr = (struct sockaddr_in *)cur->ai_addr;
        ipList << QString("%1.%2.%3.%4")
                      .arg((*addr).sin_addr.S_un.S_un_b.s_b1)
                      .arg((*addr).sin_addr.S_un.S_un_b.s_b2)
                      .arg((*addr).sin_addr.S_un.S_un_b.s_b3)
                      .arg((*addr).sin_addr.S_un.S_un_b.s_b4);
    }
    freeaddrinfo(res);
    return ipList;
}

double Round(double pSource, int pPrecision)
{
    double tmp = pow(10, pPrecision);
    return round(pSource * tmp) / tmp;
}

QString GetDir(QString pFilePath)
{
    auto pos = pFilePath.replace("\\", "/").lastIndexOf("/");
    if (pos < 0)
    {
        return "";
    }
    return pFilePath.left(pos);
}

QString SplitDir(QString &pFilePath)
{
    auto pos = pFilePath.replace("\\", "/").lastIndexOf("/");
    if (pos < 0)
    {
        return "";
    }
    auto dir = pFilePath.left(pos);
    pFilePath = pFilePath.mid(pos+1);
    return dir;
}

QString KeyValueToKeyName(int pKeyValue)
{
    QMetaEnum metaKey = QMetaEnum::fromType<Qt::Key>();
    QString keyName = metaKey.valueToKey(pKeyValue);
    return keyName.mid(4);
}

int KeyNameToKeyValue(QString pKeyName)
{
    QMetaEnum metaKey = QMetaEnum::fromType<Qt::Key>();
    QString keyName = "Key_" + pKeyName;
    return metaKey.keysToValue(keyName.toStdString().data());
}

void Wait(int pMS)
{
    QTime dieTime = QTime::currentTime().addMSecs(pMS);
    while( QTime::currentTime() < dieTime )
    {
        QApplication::processEvents(QEventLoop::AllEvents);
    }
}

QString ReadFile(QString pFile)
{
    QFile file(pFile);
    if(!file.exists())
    {
        return "";
    }
    if(file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        QString detail(file.readAll());
        file.close();
        return detail;
    }
    return "";
}

QFileInfoList GetFileList(QString pDir, QString pSuffix)
{
    QStringList filter;
    filter<<QString("*.%1").arg(pSuffix);
    return GetFileList(pDir, filter);
}

QFileInfoList GetFileList(QString pDir, QStringList pNameFilters)
{
    QDir dir(pDir);
    dir.setFilter(QDir::Files);
    QFileInfoList fileInfoList = dir.entryInfoList(pNameFilters);
    return fileInfoList;
}

}  // namespace SystemUtil
