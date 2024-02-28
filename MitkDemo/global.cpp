#include "global.h"
#include <QStorageInfo>
#include <QDirIterator>
#include <QCryptographicHash>
#include <QRegularExpression>
int AUTHORITY_USER=0;
int AUTHORITY_ADMIN=1;
QString gCustomerName = "";
user_info gCurrentUser={};
casetable_tuple gCurrentCase = {};
int gCaseId = -1;
QString gCompressPassword="compress";
QString gSpecial = "@#$%^&*";

bool IsValidId(QString pId)
{
    QRegularExpression re("^[^_@#$%^&*\\u4e00-\\u9fa5]{1,20}$");
    QRegularExpressionMatch match = re.match(pId);
    return match.hasMatch();
}

bool IsValidPassword(QString pPassword)
{
    QRegularExpression re("^(?=.*[A-Za-z])(?=.*\\d)[A-Za-z\\d]{6,}(?<![\\W_])$");
    QRegularExpressionMatch match = re.match(pPassword);
    return match.hasMatch();
}

bool IsEnoughSpace(QString pSourceDirPath, QString pTargetDirPath,quint64 pExtraSize)
{
    QDir sourceDir(pSourceDirPath);
    QDir targetDir(pTargetDirPath);

    // 获取源目录的大小
    quint64 sourceSize = 0;
    QFileInfoList fileList = sourceDir.entryInfoList(QDir::Files | QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name);
    foreach (const QFileInfo& fileInfo, fileList)
    {
        if (fileInfo.isDir())
        {
            // 子目录
            QDir childDir(fileInfo.filePath());
            sourceSize += CalculateTotalSize(childDir.path());
        }
        else
        {
            // 文件
            sourceSize += fileInfo.size();
        }
    }

    // 获取目标目录的可用大小
    QStorageInfo storage(pTargetDirPath);
    qint64 bytesAvailable = storage.bytesAvailable();

    // 比较源目录和目标目录的大小
    bool enoughSpace = (static_cast<qint64>(sourceSize+pExtraSize) < bytesAvailable);
    return enoughSpace;
}

quint64 CalculateTotalSize(const QString &pDirPath)
{
    QDirIterator it(pDirPath, QDir::Files | QDir::Dirs | QDir::NoDotAndDotDot, QDirIterator::Subdirectories);
    quint64 totalSize = 0;
    while (it.hasNext())
    {
        it.next();
        QFileInfo fileInfo = it.fileInfo();
        if (fileInfo.isFile())
        {
            totalSize += fileInfo.size();
        }
    }
    return totalSize;
}

QString encryptPassword(QString password)
{
    QString pwdmd5;
    QByteArray bb;
    bb = QCryptographicHash::hash(password.toLatin1(), QCryptographicHash::Md5);
    pwdmd5.append(bb.toHex());
    return pwdmd5;
}
