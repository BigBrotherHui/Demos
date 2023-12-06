#ifndef PLUGINEXAMPLE_GLOBAL_H
#define PLUGINEXAMPLE_GLOBAL_H

#ifdef PLUGINEXAMPLE_DLL
#define PLUGINEXAMPLE_EXPORT __declspec(dllexport)
#else
#define PLUGINEXAMPLE_EXPORT __declspec(dllimport)
#endif

#endif
