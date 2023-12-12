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
    ui->listWidget->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(ui->listWidget, &QListWidget::customContextMenuRequested, this, &Widget::on_listWidget_customContextMenuRequested);
    QList<int> sizes;
    sizes << 100 << 1000;
    ui->splitter->setSizes(sizes);
    connect(ui->pushButton_broadcast, &QPushButton::clicked, this, &Widget::on_pushButton_broadcast_clicked);
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
    for (int i = 0; i < ui->stackedWidget->count(); ++i)
    {
        ui->stackedWidget->removeWidget(ui->stackedWidget->widget(i));
    }
    ui->listWidget->clear();
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

void Widget::slot_listWidget_clicked(const QModelIndex& index)
{
    if (!index.isValid())
        return;
    QString item = ui->listWidget->item(index.row())->text();
    switchWidget(item);
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
