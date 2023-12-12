#ifndef FILE_SYSTEM_WATCHER_H
#define FILE_SYSTEM_WATCHER_H

#include <QObject>
#include <QMap>
#include <QFileSystemWatcher>
#include <QDateTime>
class FileSystemWatcher : public QObject
{
    Q_OBJECT
public:
    static FileSystemWatcher* instance();
    ~FileSystemWatcher();
    void addWatchPath(QString path);
public slots:
    void directoryUpdated(const QString& path);  // 目录更新时调用，path是监控的路径
    void fileUpdated(const QString& path);   // 文件被修改时调用，path是监控的路径
signals:
    void signal_directoryUpdated(int optype,const QStringList &names);//0添加、1删除、2覆盖
private:
    explicit FileSystemWatcher(QObject* parent = 0);
    QMap<QString,QDateTime> m_file_update_time_map;
private:
    static FileSystemWatcher* m_pInstance; // 单例
    QFileSystemWatcher* m_pSystemWatcher{nullptr},*m_pContentsWatcher{nullptr};  // QFileSystemWatcher变量
    QMap<QString, QStringList> m_currentContentsMap; // 当前每个监控的内容目录列表
    class GarbageCollector
    {
        ~GarbageCollector()
        {
            if (m_pInstance)
            {
                delete m_pInstance;
                m_pInstance = nullptr;
            }
        }
    };
    static GarbageCollector gc;
};

#endif // FILE_SYSTEM_WATCHER_H
