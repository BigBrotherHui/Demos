#ifndef CORE_GLOBAL_H
#define CORE_GLOBAL_H

#ifdef CORE_DLL
#define CORE_EXPORT __declspec(dllexport)
#else
#define CORE_EXPORT __declspec(dllimport)
#endif

#endif
