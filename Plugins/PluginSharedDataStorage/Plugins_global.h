#ifndef PluginSharedDataStorage_GLOBAL_H
#define PluginSharedDataStorage_GLOBAL_H

#ifdef PlUGINSHAREDDAtASTORAGE_DLL
#define PlUGINSHAREDDAtASTORAGE_EXPORT __declspec(dllexport)
#else
#define PlUGINSHAREDDAtASTORAGE_EXPORT __declspec(dllimport)
#endif

#endif
