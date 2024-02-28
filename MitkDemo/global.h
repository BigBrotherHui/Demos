#ifndef GLOBAL_H
#define GLOBAL_H
#pragma once
#include "QString"
#include "SqlDemo.h"
extern int AUTHORITY_ADMIN;
extern int AUTHORITY_USER;
//系统类型
extern QString gSysType;
//授权类型
extern int gLicenseType;
//序列号
extern QString gSerialNo;
//授权号
extern QString gLicense;
//过期时间
extern QString gExpireDate;
//用户名称
extern QString gCustomerName;

extern user_info gCurrentUser;
//当前已经打开的病案
extern  casetable_tuple gCurrentCase;
//当前打开（包括正在打开）的病案ID
extern int gCaseId;
extern QString gCompressPassword;
extern QString gSpecial;

extern bool IsValidId(QString pId);
extern bool IsValidPassword(QString pPassword);
extern bool IsEnoughSpace(QString pSourceDir,QString pTargetDir,quint64 pExtraSize=0);
extern quint64 CalculateTotalSize(const QString& pDirPath);
extern QString encryptPassword(QString password);
#endif // GLOBAL_H
