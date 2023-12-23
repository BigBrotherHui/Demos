#ifndef PROPERTYBROWSER_GLOBAL_H
#define PROPERTYBROWSER_GLOBAL_H

#ifdef PROPERTYBROWSER_DLL
#define PROPERTYBROWSER_EXPORT __declspec(dllexport)
#else
#define PROPERTYBROWSER_EXPORT __declspec(dllimport)
#endif

#endif
