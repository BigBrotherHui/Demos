
#include "widget.h"
#include <QApplication>
#include <QmitkRegisterClasses.h>

int main( int argc, char *argv[] )
{
	QmitkRegisterClasses();
	QApplication a(argc, argv);
	Widget w;
	w.show();
	return a.exec();
}

