//定义创建vtk智能指针对象的写法
#define vtkCreateMacro(type, obj) \
  vtkSmartPointer<type> obj = vtkSmartPointer<type>::New()

//创建单例类（需要把构造函数声明为非public，并在构造函数里DECLARE_SINGLETON(MyClass);
#define DECLARE_SINGLETON(Class) \
 Q_DISABLE_COPY(Class) \
 public: \
     static Class* Instance() \
     { \
         static QMutex mutex; \
         static QScopedPointer<Class> inst; \
         if (Q_UNLIKELY(!inst)) { \
             mutex.lock(); \
             if (!inst) inst.reset(new Class); \
             mutex.unlock(); \
         } \
         return inst.data(); \
     }

#include <qcompilerdetection.h>
#include <QMutex>
#include <qscopedpointer.h>
//声明单例类的模板方法，用法：MyClass* inst = Singleton<MyClass>::Instance();
template < class T>
class Singleton
{
public:
	static T * Instance()
	{
	    static QMutex mutex;
	    static QScopedPointer<T> inst;
	    if (Q_UNLIKELY(!inst)) {//Q_UNLIKELY告诉编译器该情况不太可能发生，使编译器优化这个代码段，与之相反的是Q_LIKELY ，是比较可能发生的情况
		     mutex.lock();
		     if (!inst) {
		          inst.reset(new T);
			 }
		     mutex.unlock();
		}
	    return inst.data();
	}
};