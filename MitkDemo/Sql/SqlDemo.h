#pragma once

#pragma execution_character_set("utf-8")
// qt
#include <QVariant>
#include <QDateTime>
#include <QSqlDatabase>
#include <QSqlQuery>
#include <QSqlError>
//
/**
 * @brief 管理员账号
 */
static const QString LOGIN_ID_ADMIN = "admin";
/**
 * @brief 医生账号
 */
static const QString LOGIN_ID_USER = "user";
/**
 * @brief 默认密码
 */
static const QString DEFAULT_PASSWORD = "123456";
enum Sex
{
    Sex_Male,
    Sex_Female
};
enum Authority
{
    Authority_User,
    Authority_Admin
};
enum Enable
{
    Enable_Able,
    Enable_Disable
};

struct casetable_tuple
{
    casetable_tuple(){}
    int case_id;
    QString case_dir;
    QDateTime case_create_date;
    QString case_operator;
    QString patient_name;
    QString patient_id;
    QString patient_gender;
    int patient_age;
    QString surgery_operator;//对应于用户id,为usertable的外键
};

struct user_info
{
    QString id;
    QString name;
    QString password;
    int sex;
    int authority;
    int enable;
};

class SqlDemo
{

public:
    static QString currentDbName;
    static bool dbExist;
    QString s_storageDir ="";

    QSqlDatabase db ;

    bool Initialize();

    bool InsertCaseTableTuple(casetable_tuple tuple);
    bool UpdateCaseTableTuple(casetable_tuple tuple);
    bool UpdateCaseTableTuple(int case_ID, QString key, QVariant value);
    bool RemoveCaseTableTuple(int case_id);
    bool AddDefaultUsers();
    bool UpdateUserPasswordById(QString pId,QString pPassword);

    int GenerateUniqueCaseID(); //return -1, if max unique cannot be Generated
    bool DisConnectDB();
    bool ConnectDB();

    QList<casetable_tuple> GetAllCaseTuple(bool *ok = nullptr);
    QList<casetable_tuple> GetCaseTupleByCaseOperator(QString pOperator,bool *ok = nullptr);
    casetable_tuple GetCaseTuple(int case_id, bool *ok = nullptr);

    //user
    bool SelectUser(QString pId,QString pPassword,user_info &pInfo);
    bool IsUserExisted(QString pId,bool &pIsExisted);
    bool AddUser(const user_info &info);
    bool GetUserByAuthority(std::vector<user_info> &pInfo,int pAuthority);
    bool GetUserByEnable(std::vector<user_info> &pInfo,Enable pEnable);
    bool RemoveUserById(QString pId);
    bool SetUserEnableById(QString pId,int pEnable=0);
    bool IsUserEnable(QString pId,int &pEnable);

//    QString SelectLoginTable(QString name,int& authority);
//    user_info SelectLoginTable(QString name);
//    bool InsertLoginTableTuple(QString name, QString password);

    bool UpdataLoginTableTuple(QString name, QString password);
    bool UpdataLoginTableTuple(QString name, QString password,bool encrypt=false);
    QString str_userpassword;
    QString SelectHash0(int name);
    QString SelectHash1(int name);
    //case
    std::vector<casetable_tuple> SearchCaseByKeyword(QString keyword);
    std::vector<casetable_tuple> SearchCaseByKeywordAndCaseOperator(QString keyword,QString pCaseOperator);
    bool RemoveCaseByCaseId(QString pCaseId);
    //settings
    QMap<QString,QString> SelectAllUserInfo();
    static QString encryptPassword(QString password);
};