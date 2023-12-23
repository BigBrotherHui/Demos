#include "globalsignalinstance.h"

GlobalSignalInstance::GlobalSignalInstance(QObject *parent) : QObject(parent) {}

GlobalSignalInstance *GlobalSignalInstance::instance() {
  static GlobalSignalInstance ins;
  return &ins;
}
