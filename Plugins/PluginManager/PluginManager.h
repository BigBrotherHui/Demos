#ifndef PLUGINMANAGER_H
#define PLUGINMANAGER_H

#include <QObject>
#include <QHash>
#include "PluginInterface.h"
#include "Plugins_global.h"
class QPluginLoader;

class PLUGINMANAGER_EXPORT PluginManager : public QObject
{
    Q_OBJECT

public:
    static PluginManager* instance();

    void loadAllPlugins();
    void loadPlugin(const QString &filepath);
    void unloadPlugin(const QString &filepath);
    void unloadAllPlugins();
    QPluginLoader* getPlugin(const QString &name);
    QVariant getPluginName(QPluginLoader *loader);
	void getAllPlugins(std::vector<QPluginLoader*> &ret);
public slots:
    void recMsgfromPlugin(PluginMetaData metadata);

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
};

#endif // PLUGINMANAGER_H
