#ifndef PLUGINMANAGER_H
#define PLUGINMANAGER_H

#include <QObject>
#include <QHash>
#include "PluginInterface.h"
#include "Core_global.h"
#include "FileSystemWatcher.h"
class QPluginLoader;
class QFileSystemWatcher;
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
    QVariant getPluginName(QPluginLoader *loader);
	void getAllPlugins(std::vector<QPluginLoader*> &ret);
public slots:
    void recMsgfromPlugin(PluginMetaData metadata);
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
};

#endif // PLUGINMANAGER_H
