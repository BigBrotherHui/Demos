//���崴��vtk����ָ������д��
#define vtkCreateMacro(type, obj) \
  vtkSmartPointer<type> obj = vtkSmartPointer<type>::New()

//���������ࣨ��Ҫ�ѹ��캯������Ϊ��public�����ڹ��캯����DECLARE_SINGLETON(MyClass);
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
//�����������ģ�巽�����÷���MyClass* inst = Singleton<MyClass>::Instance();
template < class T>
class Singleton
{
public:
	static T * Instance()
	{
	    static QMutex mutex;
	    static QScopedPointer<T> inst;
	    if (Q_UNLIKELY(!inst)) {//Q_UNLIKELY���߱������������̫���ܷ�����ʹ�������Ż��������Σ���֮�෴����Q_LIKELY ���ǱȽϿ��ܷ��������
		     mutex.lock();
		     if (!inst) {
		          inst.reset(new T);
			 }
		     mutex.unlock();
		}
	    return inst.data();
	}
};