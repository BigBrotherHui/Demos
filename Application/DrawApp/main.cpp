#include <QApplication>
#include <QDesktopWidget>
#include <QTextCodec>
#include "mainwindow.h"
int main(int argc, char *argv[]) {
  Q_INIT_RESOURCE(app);
  QApplication a(argc, argv);

  QTextCodec::setCodecForLocale(QTextCodec::codecForName("UTF-8"));
  QString fontPath_Medium = QString("%1/SourceHanSansCN-Medium.otf")
                                .arg(QCoreApplication::applicationDirPath());
  QString fontPath_Light = QString("%1/SourceHanSansCN-Light.otf")
                               .arg(QCoreApplication::applicationDirPath());
  QString fontPath_ExtraLight =
      QString("%1/SourceHanSansCN-ExtraLight.otf")
          .arg(QCoreApplication::applicationDirPath());
  QString fontPath_Bold = QString("%1/SourceHanSansCN-Bold.otf")
                              .arg(QCoreApplication::applicationDirPath());
  QString fontPath_Regular = QString("%1/SourceHanSansCN-Regular.otf")
                                 .arg(QCoreApplication::applicationDirPath());
  QString fontPath_Heavy = QString("%1/SourceHanSansCN-Heavy.otf")
                               .arg(QCoreApplication::applicationDirPath());
  QFontDatabase::addApplicationFont(fontPath_Medium);
  QFontDatabase::addApplicationFont(fontPath_Light);
  QFontDatabase::addApplicationFont(fontPath_ExtraLight);
  QFontDatabase::addApplicationFont(fontPath_Bold);
  QFontDatabase::addApplicationFont(fontPath_Regular);
  QFontDatabase::addApplicationFont(fontPath_Heavy);
  MainWindow w;
  w.showMaximized();
  return a.exec();
}
