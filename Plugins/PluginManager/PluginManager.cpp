#include "PluginManager.h"
#include <QPluginLoader>
#include <QDir>
#include <QDebug>
#include <QApplication>

PluginManager *PluginManager::m_instance=nullptr;

PluginManager::PluginManager(QObject *parent) : QObject(parent)
{

}

PluginManager::~PluginManager()
{
    unloadAllPlugins();
}

PluginManager* PluginManager::instance()
{    
	if (!m_instance)
        m_instance = new PluginManager();
	return m_instance;
}

//加载所有插件
void PluginManager::loadAllPlugins()
{
    unloadAllPlugins();
    QDir pluginsdir(qApp->applicationDirPath());
    pluginsdir.cd("plugins");

    QFileInfoList pluginsInfo = pluginsdir.entryInfoList(QDir::Files | QDir::NoDotAndDotDot);

    //加载插件
    for(QFileInfo fileinfo : pluginsInfo)
    {
        loadPlugin(fileinfo.absoluteFilePath());
    }
}

//加载其中某个插件
void PluginManager::loadPlugin(const QString &filepath)
{
    if(!QLibrary::isLibrary(filepath))
        return;
    //加载插件
    QPluginLoader *loader = new QPluginLoader(filepath);
    QString plugin_name;
    if(loader->load())
    {
        PluginInterface *plugin = dynamic_cast<PluginInterface *>(loader->instance());
        if(plugin)
        {
            plugin_name = plugin->get_name();
            m_loaders.insert(filepath, loader);
            m_names.insert(filepath,plugin_name);
            qDebug()<<"插件名称："<<plugin->get_name()<<"插件信息："<<plugin->show_text();

            connect(loader->instance(),SIGNAL(sendMsg2Manager(PluginMetaData)),this,SLOT(recMsgfromPlugin(PluginMetaData)));
        }
        else
        {
            delete loader;
            loader = nullptr;
        }
    }
    else
    {
        qDebug()<<"loadPlugin:"<<filepath<<loader->errorString();
    }
}

//卸载所有插件
void PluginManager::unloadAllPlugins()
{
    for(QString filepath : m_loaders.keys())
        unloadPlugin(filepath);
}

void PluginManager::unloadPlugin(const QString &filepath)
{
    QPluginLoader *loader = m_loaders.value(filepath);
    //卸载插件，并从内部数据结构中移除
    if(loader->unload())
    {
        m_loaders.remove(filepath);
        delete loader;
        loader = nullptr;
    }
}

//获取某个插件名称
QVariant PluginManager::getPluginName(QPluginLoader *loader)
{
    if(loader)
        return m_names.value(m_loaders.key(loader));
    else
        return "";
}

void PluginManager::getAllPlugins(std::vector<QPluginLoader*> &ret)
{
    ret.clear();
    for(auto iter=m_loaders.begin();iter!=m_loaders.end();++iter)
    {
		ret.push_back(iter.value());
    }
}

//根据名称获得插件
QPluginLoader *PluginManager::getPlugin(const QString &name)
{
    return m_loaders.value(m_names.key(name));
}

void PluginManager::recMsgfromPlugin(PluginMetaData metadata)
{
    auto loader = getPlugin(metadata.dest);//目标插件
    if(loader)
    {
        auto interface = dynamic_cast<PluginInterface*>(loader->instance());;
        if(interface)
        {
            interface->recMsgfromManager(metadata);//转发给对应插件
        }
    }
}
