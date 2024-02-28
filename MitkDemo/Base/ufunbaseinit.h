#ifndef UFUNBASEINIT_H
#define UFUNBASEINIT_H


class uFunBaseInit
{
public:
    static uFunBaseInit *getInStance();
    void f_InitBaseWidgetName();
    
private:
    explicit uFunBaseInit();
    static uFunBaseInit *self;
};

#endif // UFUNBASEINIT_H
