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
    void directoryUpdated(const QString& path);  // Ŀ¼����ʱ���ã�path�Ǽ�ص�·��
    void fileUpdated(const QString& path);   // �ļ����޸�ʱ���ã�path�Ǽ�ص�·��
signals:
    void signal_directoryUpdated(int optype,const QStringList &names);//0��ӡ�1ɾ����2����
private:
    explicit FileSystemWatcher(QObject* parent = 0);
    QMap<QString,QDateTime> m_file_update_time_map;
private:
    static FileSystemWatcher* m_pInstance; // ����
    QFileSystemWatcher* m_pSystemWatcher{nullptr},*m_pContentsWatcher{nullptr};  // QFileSystemWatcher����
    QMap<QString, QStringList> m_currentContentsMap; // ��ǰÿ����ص�����Ŀ¼�б�
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
