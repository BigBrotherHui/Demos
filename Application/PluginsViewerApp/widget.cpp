#include "widget.h"
#include "ui_widget.h"
#include <QDir>
#include <QDebug>
#include <QPluginLoader>
#include "PluginInterface.h"
#include "PluginManager.h"
#include <QModelIndex>
Widget::Widget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Widget)
{
    ui->setupUi(this);
    QList<int> sizes;
    sizes << 100 << 1000;
    ui->splitter->setSizes(sizes);
    connect(ui->pushButton_broadcast, &QPushButton::clicked, this, [&]
        {
            for(int i=0;i<ui->stackedWidget->count();++i)
            {
                ui->stackedWidget->removeWidget(ui->stackedWidget->widget(i));
            }
            ui->listWidget->clear();
            PluginManager* pm = PluginManager::instance();
            std::vector<QPluginLoader*> plugins;
    		pm->getAllPlugins(plugins);
            for(int i=0;i<plugins.size();++i)
            {
                QPluginLoader* loader = plugins.at(i);
                if (loader)
                {
                    PluginUIInterface* plugin = dynamic_cast<PluginUIInterface*>(loader->instance());
                    if (plugin && plugin->createWidget())
                    {
                        ui->listWidget->addItem(plugin->get_name());
                        ui->stackedWidget->insertWidget(0, plugin->createWidget());
                    }
                }
                else
                    qDebug() << "插件不存在";
            }
        });
    connect(ui->listWidget, &QListWidget::clicked, this, [&](const QModelIndex& index)
        {
            if (!index.isValid())
                return;
            QString item = ui->listWidget->item(index.row())->text();
            ui->stackedWidget->setCurrentWidget(dynamic_cast<PluginUIInterface*>(PluginManager::instance()->getPlugin(item)->instance())->createWidget());
        });
}

Widget::~Widget()
{
    delete ui;
}

//加载所有插件
void Widget::on_pushButton_clicked()
{
    PluginManager * pm = PluginManager::instance();
    pm->loadAllPlugins();
}

//卸载所有插件
void Widget::on_pushButton_3_clicked()
{
    PluginManager * pm = PluginManager::instance();
    pm->unloadAllPlugins();
}

//插件01发消息
void Widget::on_pushButton_2_clicked()
{
    PluginManager * pm = PluginManager::instance();
    auto loader = pm->getPlugin("Plugin01");
    if(loader)
    {
        PluginInterface *plugin = dynamic_cast<PluginInterface *>(loader->instance());
        PluginMetaData m;
        m.dest = "Plugin02";
        m.from = "Plugin01";
        m.msg = "插件1给插件2发的消息";
        plugin->sendMsg2Manager(m);
    }
    else
        qDebug()<<"插件不存在";
}

void Widget::on_pushButton_4_clicked()
{
    PluginManager * pm = PluginManager::instance();
    auto loader = pm->getPlugin("Plugin02");
    if(loader)
    {
        PluginInterface *plugin = dynamic_cast<PluginInterface *>(loader->instance());
        PluginMetaData m;
        m.dest = "Plugin01";
        m.from = "Plugin02";
        m.msg = "插件2给插件1发的消息";
        emit plugin->sendMsg2Manager(m);
    }
}

void Widget::on_pushButton_20_clicked()
{
    PluginManager * pm = PluginManager::instance();
    auto loader = pm->getPlugin("Plugin_widget");
    if(loader)
    {
        PluginInterface *plugin = dynamic_cast<PluginInterface *>(loader->instance());
        PluginMetaData m;
        m.dest = "Plugin_widget";
        m.from = "";
        m.msg = "show";
        emit plugin->sendMsg2Manager(m);
    }
}
