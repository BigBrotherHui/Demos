
#include "widget.h"
//#include "registerfineadjustdialog.h"
#include <QApplication>
#include <QmitkRegisterClasses.h>

int main( int argc, char *argv[] )
{
	QmitkRegisterClasses();
	QApplication a(argc, argv);
	//RegisterFineAdjustDialog w;
	Widget w;
	w.show();
	return a.exec();
}

