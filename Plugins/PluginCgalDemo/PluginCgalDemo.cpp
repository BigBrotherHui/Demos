#include <QDir>
#include <QDebug>
#include <QPluginLoader>
#include "PluginInterface.h"
#include "PluginCgalDemo.h"
#include "Widget.h"
int PluginCgalDemo::typeId = qRegisterMetaType<PluginCgalDemo*>("PluginCgalDemo");
PluginCgalDemo::PluginCgalDemo()
{
}

PluginCgalDemo::~PluginCgalDemo()
{
}

QWidget* PluginCgalDemo::createWidget()
{
    if (!m_widget)
        m_widget = new Widget(nullptr);
    return m_widget;
}