#ifndef PLUGINMANAGER_GLOBAL_H
#define PLUGINMANAGER_GLOBAL_H

#ifdef PLUGINMANAGER_DLL
#define PLUGINMANAGER_EXPORT __declspec(dllexport)
#else
#define PLUGINMANAGER_EXPORT __declspec(dllimport)
#endif

#endif
