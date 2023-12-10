#ifndef PLUGINCGALDEMO_GLOBAL_H
#define PLUGINCGALDEMO_GLOBAL_H

#ifdef PLUGINCGALDEMO_DLL
#define PLUGINCGALDEMO_EXPORT __declspec(dllexport)
#else
#define PLUGINCGALDEMO_EXPORT __declspec(dllimport)
#endif

#endif
