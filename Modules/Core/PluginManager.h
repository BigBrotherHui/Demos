#ifndef PLUGINMANAGER_H
#define PLUGINMANAGER_H

#include <QObject>
#include <QHash>
#include "PluginInterface.h"
#include "Core_global.h"
#include "FileSystemWatcher.h"
#include <queue>
#include <QMutex>
#include <QWaitCondition>
#include <QThread>
class QPluginLoader;
class QFileSystemWatcher;
struct compare //重写仿函数
{
    bool operator() (const PluginMetaData& a, const PluginMetaData& b)
    {
        return a.priority < b.priority; //大顶堆
    }
};
class SendThreadObject : public QObject
{
    Q_OBJECT
public:
    SendThreadObject(QMutex* mutex, QWaitCondition* condition, std::priority_queue<PluginMetaData, std::vector<PluginMetaData>, compare>* queue);
public slots:
    void slot_DealMetadata();
private:
    QMutex* m_mutex;
    QWaitCondition* m_condition;
    std::priority_queue<PluginMetaData, std::vector<PluginMetaData>, compare>* m_queue;
};
class CORE_EXPORT PluginManager : public QObject
{
    Q_OBJECT

public:
    static PluginManager* instance();
    QObject* tryConstructObject(QString pageName);
    void loadAllPlugins();
    bool loadPlugin(const QString &filepath);
    bool unloadPlugin(const QString &filepath);
    bool unloadPlugin(QPluginLoader* loader);
    void unloadAllPlugins();
    QPluginLoader* getPlugin(const QString &name);
    QString getPluginName(QPluginLoader *loader);
	void getAllPlugins(std::vector<QPluginLoader*> &ret);
public slots:
    void recMsgfromPlugin(PluginMetaData metadata);
protected:
    void sendMsg(const PluginMetaData& data);
protected slots:
    void slot_directoryUpdated(int optype,const QStringList &names);
private:
    explicit PluginManager(QObject *parent = nullptr);
    ~PluginManager();
    static PluginManager* m_instance;
    QHash<QString, QPluginLoader *> m_loaders; //插件路径--QPluginLoader实例
    QHash<QString, QString> m_names; //插件路径--插件名称
    class GarbageCollector
    {
        ~GarbageCollector()
        {
            if (m_instance)
            {
                delete m_instance;
                m_instance = nullptr;
            }
        }
    };
    static GarbageCollector gc;
    FileSystemWatcher *m_file_system_watcher;
    std::priority_queue<PluginMetaData, std::vector<PluginMetaData>, compare> m_waiting_metadatas;
    QMutex m_mutex;
    QWaitCondition m_condition;
    friend class SendThreadObject;
    SendThreadObject *m_sendObject;
    QThread* m_sendThread;
};

#endif // PLUGINMANAGER_H
