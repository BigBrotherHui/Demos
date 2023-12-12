#include <QDir>
#include <QFileInfo>
#include <qDebug>
#include "FileSystemWatcher.h"

FileSystemWatcher* FileSystemWatcher::m_pInstance = NULL;

FileSystemWatcher::FileSystemWatcher(QObject* parent)
    : QObject(parent)
{

}

FileSystemWatcher* FileSystemWatcher::instance()
{
    if (m_pInstance == NULL)
    {
        m_pInstance = new FileSystemWatcher();
        m_pInstance->m_pSystemWatcher = new QFileSystemWatcher();
        m_pInstance->m_pContentsWatcher = new QFileSystemWatcher();
        // ����QFileSystemWatcher��directoryChanged��fileChanged�źŵ���Ӧ�Ĳ�
        connect(m_pInstance->m_pSystemWatcher, SIGNAL(directoryChanged(QString)), m_pInstance, SLOT(directoryUpdated(QString)));
        connect(m_pInstance->m_pContentsWatcher, SIGNAL(fileChanged(QString)), m_pInstance, SLOT(fileUpdated(QString)));
    }
    return m_pInstance;
}

FileSystemWatcher::~FileSystemWatcher()
{
    if(m_pSystemWatcher)
    {
        m_pSystemWatcher->deleteLater();
        m_pSystemWatcher = nullptr;
    }
    if (m_pContentsWatcher)
    {
        m_pContentsWatcher->deleteLater();
        m_pContentsWatcher = nullptr;
    }
}

// ����ļ���Ŀ¼
void FileSystemWatcher::addWatchPath(QString path)
{
    qDebug() << QString("Add to watch: %1").arg(path);
    // ��Ӽ��·��
    m_pSystemWatcher->addPath(path);
    
    // ������·����һ��Ŀ¼�����浱ǰ�����б�
    QFileInfo file(path);
    if (file.isDir())
    {
        const QDir dirw(path);
        for (auto p : m_currentContentsMap[path])
        {
            m_pContentsWatcher->removePath(path + "/" + p);
        }
        m_currentContentsMap[path] = dirw.entryList(QDir::NoDotAndDotDot | QDir::AllDirs | QDir::Files, QDir::DirsFirst);
        for(auto p : m_currentContentsMap[path])
        {
            m_pContentsWatcher->addPath(path + "/"+p);
        }
    }
}

// ֻҪ�κμ�ص�Ŀ¼���£���ӡ�ɾ���������������ͻ���á�
void FileSystemWatcher::directoryUpdated(const QString& path)
{
    // �Ƚ����µ����ݺͱ���������ҳ�����(�仯)
    QStringList currEntryList = m_currentContentsMap[path];
    const QDir dir(path);

    QStringList newEntryList = dir.entryList(QDir::NoDotAndDotDot | QDir::AllDirs | QDir::Files, QDir::DirsFirst);

    QSet<QString> newDirSet = QSet<QString>::fromList(newEntryList);
    QSet<QString> currentDirSet = QSet<QString>::fromList(currEntryList);

    // ������ļ�
    QSet<QString> newFiles = newDirSet - currentDirSet;
    QStringList newFile = newFiles.toList();
    
    // �ļ��ѱ��Ƴ�
    QSet<QString> deletedFiles = currentDirSet - newDirSet;
    QStringList deleteFile = deletedFiles.toList();

    // ���µ�ǰ����
    m_currentContentsMap[path] = newEntryList;
    for (auto p : m_currentContentsMap[path])
    {
        m_pContentsWatcher->removePath(path + "/" + p);
    }
    for (auto p : m_currentContentsMap[path])
    {
        m_pContentsWatcher->addPath(path + "/" + p);
    }
    if (!newFile.isEmpty() && !deleteFile.isEmpty())
    {
        // �ļ�/Ŀ¼������
        if ((newFile.count() == 1) && (deleteFile.count() == 1))
        {
            qDebug() << QString("File Renamed from %1 to %2").arg(deleteFile.first()).arg(newFile.first());
        }
    }
    else
    {
        // ������ļ�/Ŀ¼��Dir
        if (!newFile.isEmpty())
        {
            qDebug() << "New Files/Dirs added: " << newFile;

            foreach(QString file, newFile)
            {
                // �������ÿ�����ļ�....
            }
            emit signal_directoryUpdated(0, newFile);
        }

        // ��Dir��ɾ���ļ�/Ŀ¼
        if (!deleteFile.isEmpty())
        {
            qDebug() << "Files/Dirs deleted: " << deleteFile;

            foreach(QString file, deleteFile)
            {
                // �������ÿ����ɾ�����ļ�....
            }
            emit signal_directoryUpdated(1, deleteFile);
        }

        //����
        if(newFile.isEmpty() && deleteFile.isEmpty())
        {
            //emit signal_directoryUpdated(2, QStringList());
        }
    }
}

// �ļ��޸�ʱ����
void FileSystemWatcher::fileUpdated(const QString& path)
{
    QFileInfo file(path);
    QString strPath = file.absolutePath();
    QString strName = file.fileName();
    if(!m_file_update_time_map.contains(path))
    {
		m_file_update_time_map[path]= QDateTime::currentDateTime();
    }
	else if (m_file_update_time_map[path].msecsTo(QDateTime::currentDateTime()) < 500)
    {
		return;
    }
    m_file_update_time_map[path] = QDateTime::currentDateTime();
    emit signal_directoryUpdated(2, QStringList() << strPath + "/" + strName);
}
