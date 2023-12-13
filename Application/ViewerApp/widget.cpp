#include "widget.h"
#include "ui_widget.h"
#include <QDir>
#include <QDebug>
#include <QPluginLoader>
#include "PluginInterface.h"
#include "PluginManager.h"
#include <QModelIndex>
#include <QJsonObject>
#include <QJsonDocument>
#include <QJsonArray>
Widget::Widget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Widget)
{
    ui->setupUi(this);
    ui->listWidget->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(ui->listWidget, &QListWidget::customContextMenuRequested, this, &Widget::on_listWidget_customContextMenuRequested);
    QList<int> sizes;
    sizes << 100 << 1000;
    ui->splitter->setSizes(sizes);
    connect(ui->listWidget, &QListWidget::clicked, this, &Widget::slot_listWidget_clicked);
    m_popMenu = new QMenu(this);
    QAction* unloadAction = new QAction(tr("unload"), this);
    connect(unloadAction, SIGNAL(triggered()), this, SLOT(slot_unloadAction_triggered()));
    m_popMenu->addAction(unloadAction);
}

Widget::~Widget()
{
    delete ui;
}

void Widget::switchWidget(QString pluginName) {
    auto instance = dynamic_cast<PluginUIInterface*>(PluginManager::instance()->getPlugin(pluginName)->instance());
    if(instance)
		ui->stackedWidget->setCurrentWidget(instance->createWidget());
}

//加载所有插件
void Widget::on_pushButton_clicked()
{
    PluginManager * pm = PluginManager::instance();
    pm->loadAllPlugins();
}

void Widget::on_listWidget_customContextMenuRequested(const QPoint& pos)
{
    static QDateTime dt = QDateTime::currentDateTime();
    static bool first{1};
    if (!first && dt.msecsTo(QDateTime::currentDateTime()) < 500)
        return;
    first = false;
    dt = QDateTime::currentDateTime();
    QListWidgetItem* curItem = ui->listWidget->itemAt(pos);
    if (curItem == nullptr)
        return;
    m_popMenu->exec(QCursor::pos());
}

//卸载所有插件
void Widget::on_pushButton_3_clicked()
{
    PluginManager * pm = PluginManager::instance();
    pm->unloadAllPlugins();
    ui->listWidget->clear();
}

void Widget::on_pushButton_broadcast_clicked()
{
    bool usePageConfig{ false };
    for (int i = 0; i < ui->stackedWidget->count(); ++i)
    {
        ui->stackedWidget->removeWidget(ui->stackedWidget->widget(i));
    }
    ui->listWidget->clear();
    QString pageConfigPath = qApp->applicationDirPath() + "/" + "pageConfig.json";
    QFile file(pageConfigPath);
    if(!file.open(QIODevice::ReadWrite))
    {
        qDebug() << "pageConfig.json open failed";
        return;
    }
    QByteArray allData = file.readAll();
    QJsonParseError json_error;
    QJsonDocument jsonDoc(QJsonDocument::fromJson(allData, &json_error));
    qDebug() << "json error no:" << json_error.error;
    if (QJsonParseError::IllegalValue == json_error.error)
    {
        QJsonObject rootObj = jsonDoc.object();
        if (!rootObj.contains("PageConfig"))
        {
            QJsonArray jsonArray;
            QJsonObject page1;
            page1.insert("PageIndex", 0);
            page1.insert("PluginClassName", "PluginPolyDataToImage");
            jsonArray.append(page1);
            QJsonObject jsonObject;
            jsonObject.insert("PageConfig", jsonArray);
            jsonDoc.setObject(jsonObject);
            file.write(jsonDoc.toJson());
        }
    }
    if (json_error.error == QJsonParseError::NoError || QJsonParseError::IllegalValue == json_error.error)
    {
        QJsonObject rootObj = jsonDoc.object();
        if (rootObj.contains("PageConfig"))
        {
            QJsonValue value = rootObj.value("PageConfig");
            if (value.isArray())
            {
                QJsonArray array = value.toArray();
                size_t size = array.size();
                for (int i = 0; i < size; ++i)
                {
                    int pageIndex = array.at(i).toObject().value("PageIndex").toInt();
                    QString pluginClassName= array.at(i).toObject().value("PluginClassName").toString();
                    QObject *obj=PluginManager::instance()->tryConstructObject(pluginClassName);
                    if (!obj)
                    {
                        qDebug() << "page config error";
                        return;
                    }
                    usePageConfig = 1;
                    PluginUIInterface *interface=dynamic_cast<PluginUIInterface*>(obj);
                    if(interface)
                    {
                        ui->listWidget->addItem(interface->get_name());
                        ui->stackedWidget->insertWidget(pageIndex, interface->createWidget());
                    }
                }
            }
        }
        else
        {
            QJsonArray jsonArray;
            QJsonObject page1;
            page1.insert("page0", "PluginPolyDataToImage");
            jsonArray.append(page1);
            QJsonObject jsonObject;
            jsonObject.insert("PageConfig", jsonArray);
            jsonDoc.setObject(jsonObject);
            file.write(jsonDoc.toJson());
        }
    }
    file.close();
    if(!usePageConfig)
    {
        PluginManager* pm = PluginManager::instance();
        std::vector<QPluginLoader*> plugins;
        pm->getAllPlugins(plugins);
        for (int i = 0; i < plugins.size(); ++i)
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
    }
}

void Widget::slot_listWidget_clicked(const QModelIndex& index)
{
    if (!index.isValid())
        return;
    //QString item = ui->listWidget->item(index.row())->text();
    //switchWidget(item);
    ui->stackedWidget->setCurrentIndex(index.row());
}

void Widget::slot_unloadAction_triggered()
{
    if(PluginManager::instance()->unloadPlugin(PluginManager::instance()->getPlugin(ui->listWidget->currentItem()->text())))
    {
        ui->listWidget->takeItem(ui->listWidget->row(ui->listWidget->currentItem()));
    }else
    {
        qDebug() << "unload " << ui->listWidget->currentItem()->text() << " failed";
    }
}
