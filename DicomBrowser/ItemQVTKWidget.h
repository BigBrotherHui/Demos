#pragma once
#include<QListWidget.h>
#include <QVTKOpenGLNativeWidget.h>
namespace Ui {
	class ItemQVTKWidget;
}

class ItemQVTKWidget:public QObject , public QListWidgetItem
{
	Q_OBJECT

public:
	explicit ItemQVTKWidget(QVTKOpenGLNativeWidget*parent = 0);
	~ItemQVTKWidget();

	//…Ë÷√ ˝æ›
	void SetData(const QString& qstrFileName, int iFileSize, const QString& qstrPic);

private:

};

