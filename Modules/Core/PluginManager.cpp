#include "PluginManager.h"
#include <QPluginLoader>
#include <QDir>
#include <QDebug>
#include <QApplication>
#include <QThread>
#include <QFileSystemWatcher>
#include "UsefulFunctions.h"
PluginManager *PluginManager::m_instance=nullptr;

PluginManager::PluginManager(QObject *parent) : QObject(parent)
{
    m_file_system_watcher = FileSystemWatcher::instance();
    m_file_system_watcher->addWatchPath(qApp->applicationDirPath() + "/plugins");
    connect(m_file_system_watcher, &FileSystemWatcher::signal_directoryUpdated, this, &PluginManager::slot_directoryUpdated);
    m_sendObject = new SendThreadObject(&m_mutex,&m_condition,&m_waiting_metadatas);
    m_sendThread = new QThread(this);
    m_sendObject->moveToThread(m_sendThread);
    connect(m_sendThread, &QThread::started, m_sendObject, &SendThreadObject::slot_DealMetadata);
    m_sendThread->start();
}

PluginManager::~PluginManager()
{
    unloadAllPlugins();
    m_sendThread->quit();
    m_sendThread->wait();
}

SendThreadObject::SendThreadObject(QMutex* mutex, QWaitCondition* condition, std::priority_queue<PluginMetaData, std::vector<PluginMetaData>, compare>* queue)
{
    m_mutex = mutex;
    m_condition = condition;
    m_queue = queue;
}

void SendThreadObject::slot_DealMetadata()
{
    while(1)
    {
        if (m_queue->size() == 0)
        {
            m_mutex->lock();
            m_condition->wait(m_mutex);   /*生产者消费者模型, 主要是要理解wait的机制, wait时, 会阻塞在这里, 且同时也会把锁打开(让生产者可以向queue中插入数据).
            生产者向queue中插入数据后, 需要调用condition.wakeall()或wakeone(), 当condition收到wake信号时, 会把锁重新锁上(暂时不让生产者向queue中插数据), 并开始执行此处
            后面的代码.*/
            if(m_queue->size() > 0)
            {
                PluginMetaData meta=m_queue->top();
                PluginManager::instance()->sendMsg(meta);
                m_queue->pop();
            }
            m_mutex->unlock();
        }
    }
}

PluginManager* PluginManager::instance()
{    
	if (!m_instance)
        m_instance = new PluginManager();
	return m_instance;  
}

QObject* PluginManager::tryConstructObject(QString pageName)
{
    qDebug() << "tryConstructObject" << pageName;
    int type = QMetaType::type((pageName).toStdString().c_str());
    const QMetaObject* metaObj = QMetaType::metaObjectForType(type);
    if (metaObj)
        return metaObj->newInstance();
    else
        return nullptr;
}

//加载所有插件
void PluginManager::loadAllPlugins()
{
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
bool PluginManager::loadPlugin(const QString &filepath)
{
    if(!QLibrary::isLibrary(filepath))
        return false;
    if(m_names.contains(filepath))
    {
        qDebug() << "plugin " << filepath << " existed";
        return false;
    }
    //加载插件
    QPluginLoader *loader = new QPluginLoader(filepath);
    loader->setLoadHints(0);//只有加上这句话才可以在unload之后真正卸载该dll（即可以删除它）
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
            qDebug() << filepath << " loaded";
            return true;
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
    return false;
}

//卸载所有插件
void PluginManager::unloadAllPlugins()
{
    for(QString filepath : m_loaders.keys())
        unloadPlugin(filepath);
}

bool PluginManager::unloadPlugin(const QString &filepath)
{
    if (!QFile::exists(filepath))
        return false;
    if (!m_loaders.contains(filepath))
        return false; 
    QPluginLoader *loader = m_loaders.value(filepath);
    //卸载插件，并从内部数据结构中移除
    if(loader->unload())
    {
        m_loaders.remove(filepath);
        m_names.remove(filepath);
        delete loader;
        loader = nullptr;
        qDebug() << filepath << " unloaded";
        return true;
    }else
    {
        qDebug() << "unload failed:"<<loader->errorString();
        return false;
    }
}

bool PluginManager::unloadPlugin(QPluginLoader* loader)
{
    return unloadPlugin(m_names.key(getPluginName(loader)));
}

//获取某个插件名称
QString PluginManager::getPluginName(QPluginLoader *loader)
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

void PluginManager::sendMsg(const PluginMetaData& data)
{
    if (data.dest.isEmpty())
    {
        qDebug() << "sendMsg error:dest is empty";
        return;
    }
    auto loader = getPlugin(data.dest);
    if (!loader)
        return;
    auto interface = dynamic_cast<PluginInterface*>(loader->instance());;
    if (!interface)
        return;
    interface->recMsgfromManager(data);
}

void PluginManager::recMsgfromPlugin(PluginMetaData metadata)
{
    if(metadata.blocking)
    {
        if (!metadata.dest.isEmpty())
        {
            sendMsg(metadata);
        }
        else
        {
            for (auto loader : m_loaders)
            {
                QString name = getPluginName(loader);
                metadata.dest = name;
                sendMsg(metadata);
            }
        }
    }else
    {
        if (!metadata.dest.isEmpty())
        {
            m_mutex.lock();
            m_waiting_metadatas.push(metadata);
            m_condition.wakeAll();
            m_mutex.unlock();
        }
        else
        {
            for (auto loader : m_loaders)
            {
                QString name = getPluginName(loader);
                metadata.dest = name;
                m_mutex.lock();
                m_waiting_metadatas.push(metadata);
                m_condition.wakeAll();
                m_mutex.unlock();
            }
        }
    }
}

void PluginManager::slot_directoryUpdated(int optype, const QStringList& names)
{
	if(optype==0)//添加
	{
        for (auto name : names)
        {
            loadPlugin(name);
        }
	}else if(optype==1)//删除
	{
        for (auto name : names)
        {
            unloadPlugin(name);
        }
	}else if(optype==2)//修改
	{
        for(auto name : names)
        {
            loadPlugin(name);
        }
	}
}
