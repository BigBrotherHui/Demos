#pragma execution_character_set("utf-8")
#include "SqlDemo.h"
// qt

#include <QDir>
#include <QFile>
#include <QDebug>
#include <QDateTime>
#include <QCryptographicHash>
#include <QSettings>
#include <QCoreApplication>
#include <QFileSystemWatcher>
#include <QMessageBox>
//
#include "ustatus.h"
#include "StrSql.h"
#include <mitkLogMacros.h>

SqlDemo *mSql;

void f_UpdateCaseTableTupleFunction(int case_ID, QString key, QVariant value)
{
    mSql->UpdateCaseTableTuple(case_ID, key, value);
}

QString SqlDemo::currentDbName;
bool SqlDemo::dbExist=true;

bool SqlDemo::Initialize()
{
    s_storageDir = QCoreApplication::applicationDirPath();
    uStatus::pFun_UpdateHash = &(f_UpdateCaseTableTupleFunction);
    if (!ConnectDB())
    {
        qDebug() << "db connect fail";
        return false;
    }

    if (!db.tables().contains(QLatin1String("usertable")))
    {
        QSqlQuery sql_query;
        sql_query.prepare(SQL_CreateLoginTable);
        if (!sql_query.exec())
        {
            qDebug() << "Error : create database failed" << sql_query.lastError().text();
            return false;
        }
        AddDefaultUsers();
    }

    if (!db.tables().contains( QLatin1String("casetable")))
    {
        QSqlQuery sql_query;
        sql_query.prepare(SQL_CreateCaseTable);

        if(!sql_query.exec())
        {
            qDebug() << "Error : create database failed" << sql_query.lastError().text();
            return false;
        }
    }
  
    return true;
}


bool SqlDemo::DisConnectDB()
{
    QString name;
    {
        name = QSqlDatabase::database().connectionName();
    }//
    QSqlDatabase::removeDatabase(name);

    return true;
}

bool SqlDemo::ConnectDB()
{
    uStatus::mMainDataDir = QCoreApplication::applicationDirPath() + "/data";
    uStatus::mMainPath = QCoreApplication::applicationDirPath();
    s_storageDir=uStatus::mMainDataDir+"/db";
    //Check if storage directory exist, create one if not exist
    QDir dir_storage(s_storageDir);
    if (!dir_storage.exists()) dir_storage.mkpath(s_storageDir);
    if (!dir_storage.exists()) return EXIT_FAILURE;

    currentDbName = s_storageDir + "/MitkDemo.db";
    if(!QFile(currentDbName).exists())
    {
        if(QMessageBox::No ==QMessageBox::question(nullptr,"提示","未找到数据库文件，是否重新创建数据库？",QMessageBox::Yes|QMessageBox::No))
        {
            return false;
        }
    }
    db = QSqlDatabase::addDatabase("QSQLITE");
    db.setHostName("local");
    db.setDatabaseName(currentDbName);
    
    bool ok = db.open();
    if (!ok)
    {
        QSqlError error = db.lastError();
        QMessageBox::warning(nullptr,"提示","打开数据库失败，可能数据库已损坏，请联系售后服务人员！\n系统将退出运行！");
        return false;
    }
    return true;
}

bool SqlDemo::AddDefaultUsers()
{
    QSqlQuery sql_query;
    sql_query.prepare(SQL_InsertUser);
    sql_query.addBindValue(LOGIN_ID_ADMIN);
    sql_query.addBindValue("admin");
    sql_query.addBindValue(encryptPassword(DEFAULT_PASSWORD));
    sql_query.addBindValue(Sex_Male);
    sql_query.addBindValue(Authority_Admin);
    sql_query.addBindValue(Enable_Able);
    if (!sql_query.exec())
    {
        qDebug() << "默认管理员添加失败："<<sql_query.lastError().text();
        return 0;
    }
    sql_query.prepare(SQL_InsertUser);
    sql_query.addBindValue(LOGIN_ID_USER);
    sql_query.addBindValue("doctor");
    sql_query.addBindValue(encryptPassword(DEFAULT_PASSWORD));
    sql_query.addBindValue(Sex_Male);
    sql_query.addBindValue(Authority_User);
    sql_query.addBindValue(Enable_Able);
    if (!sql_query.exec())
    {
        qDebug() << "默认用户添加失败："<<sql_query.lastError().text();
        return 0;
    }
    return 1;
}

bool SqlDemo::UpdateUserPasswordById(QString pId,QString pPassword)
{
    QSqlQuery q;
    q.prepare(SQL_UpdateUserPasswordById);
    q.addBindValue(pPassword);
    q.addBindValue(pId);
    if(!q.exec())
        return false;
    return true;
}

bool SqlDemo::InsertCaseTableTuple(casetable_tuple tuple)
{
    QString insert_sql = "insert into casetable values (?, ?, ?, ?, ?, ?, ?, ?,?)";
    QSqlQuery sql_query;
    sql_query.prepare(insert_sql);
    sql_query.addBindValue(tuple.case_id);
    sql_query.addBindValue(tuple.case_dir);
    sql_query.addBindValue(tuple.case_create_date);
    sql_query.addBindValue(tuple.case_operator);
    sql_query.addBindValue(tuple.patient_name);
    sql_query.addBindValue(tuple.patient_id);
    sql_query.addBindValue(tuple.patient_gender);
    sql_query.addBindValue(tuple.patient_age);
    sql_query.addBindValue(tuple.surgery_operator);
    if(!sql_query.exec())
    {
        qDebug() << sql_query.lastError().text();
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

bool SqlDemo::RemoveCaseTableTuple(int case_id)
{
    QString delete_sql = "delete from casetable where case_id = ?";
    QSqlQuery sql_query;
    sql_query.prepare(delete_sql);
    sql_query.addBindValue(case_id);

    if(!sql_query.exec())
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

bool SqlDemo::UpdateCaseTableTuple(casetable_tuple tuple)
{
    //update row
    QString age = QString::number(tuple.patient_age);
    QString id = QString::number(tuple.case_id);

    QString update_sql="update casetable set patient_name=?,patient_id=?,patient_gender=?,patient_age=?,surgery_operator=?,"
                         "case_operator=? where case_id =?";
    QSqlQuery sql_query;
    sql_query.prepare(update_sql);
    sql_query.addBindValue(tuple.patient_name);
    sql_query.addBindValue(tuple.patient_id);
    sql_query.addBindValue(tuple.patient_gender);
    sql_query.addBindValue(age);
    sql_query.addBindValue(tuple.surgery_operator);
    sql_query.addBindValue(tuple.case_operator);
    sql_query.addBindValue(id);
    if(!sql_query.exec())
    {
        qDebug() << sql_query.lastError();
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

bool SqlDemo::UpdateCaseTableTuple(int case_ID, QString key, QVariant value)
{
    //update row
    QString update_sql;
    update_sql = "update casetable set {key} = :val where case_id = :c_id";
    update_sql = update_sql.replace("{key}",key);

    QSqlQuery sql_query;
    sql_query.prepare(update_sql);
    sql_query.bindValue(":val" , value);
    sql_query.bindValue(":c_id", case_ID);

    if (!sql_query.exec())
    {
        return EXIT_FAILURE;
    }

    return false;
}


int SqlDemo::GenerateUniqueCaseID()
{
    QList<int> list_id;

    QString select_sql = "select case_id from casetable";
    QSqlQuery sql_query;
    sql_query.prepare(select_sql);

    if(!sql_query.exec())
    {
        qDebug()<<sql_query.lastError();
        return -1;
    }

    int max_id = 0;
    while(sql_query.next())
    {
        if (sql_query.value(0).toInt() > max_id)
            max_id = sql_query.value(0).toInt();
    }

    return max_id + 1;
}

QList<casetable_tuple> SqlDemo::GetAllCaseTuple(bool *ok)
{
    QList<casetable_tuple> list_tuple;

    QString select_sql = "select * from casetable order by case_id DESC";          //Deafult sort by case_id in descending order
    QSqlQuery sql_query;
    sql_query.prepare(select_sql);

    if(!sql_query.exec())
    {
        qDebug()<<sql_query.lastError();
        if (ok !=nullptr)
            *ok = false;

        return list_tuple;
    }

    while(sql_query.next())
    {
        int i=0;
        casetable_tuple tp;
        tp.case_id           = sql_query.value(i++).toInt();
        tp.case_dir          = sql_query.value(i++).toString();
        tp.case_create_date  = sql_query.value(i++).toDateTime();
        tp.case_operator=sql_query.value(i++).toString();
        tp.patient_name      = sql_query.value(i++).toString();
        tp.patient_id        = sql_query.value(i++).toString();
        tp.patient_gender    = sql_query.value(i++).toString();
        tp.patient_age       = sql_query.value(i++).toInt();
        tp.surgery_operator  = sql_query.value(i++).toString();
        list_tuple.append(tp);
    }

    return list_tuple;
}

QList<casetable_tuple> SqlDemo::GetCaseTupleByCaseOperator(QString pOperator, bool *ok)
{
    QList<casetable_tuple> list_tuple;
    QString select_sql = QString("select * from casetable where surgery_operator=? order by case_id DESC");  
    QSqlQuery sql_query;
    sql_query.prepare(select_sql);
    sql_query.addBindValue(pOperator);
    if(!sql_query.exec())
    {
        qDebug()<<sql_query.lastError();
        if (ok !=nullptr)
            *ok = false;
        return list_tuple;
    }
    while(sql_query.next())
    {
        int i=0;
        casetable_tuple tp;
        tp.case_id           = sql_query.value(i++).toInt();
        tp.case_dir          = sql_query.value(i++).toString();
        tp.case_create_date  = sql_query.value(i++).toDateTime();
        tp.case_operator=sql_query.value(i++).toString();
        tp.patient_name      = sql_query.value(i++).toString();
        tp.patient_id        = sql_query.value(i++).toString();
        tp.patient_gender    = sql_query.value(i++).toString();
        tp.patient_age       = sql_query.value(i++).toInt();
        tp.surgery_operator  = sql_query.value(i++).toString();
        list_tuple.append(tp);
    }
    return list_tuple;
}

casetable_tuple SqlDemo::GetCaseTuple(int case_id, bool *ok)
{
    casetable_tuple selectedTuple;
    QString select_sql = "select * from casetable where case_id = ?";
    QSqlQuery sql_query;
    sql_query.prepare(select_sql);
    sql_query.addBindValue(case_id);

    //return if query fail
    if(!sql_query.exec())
    {
        qDebug()<<sql_query.lastError();
        if (ok !=nullptr)
            *ok = false;
        return selectedTuple;
    }

    if(!sql_query.next())
    {
        if (ok !=nullptr)
            *ok = false;
        return selectedTuple;
    }

    selectedTuple.case_id           = sql_query.value(0).toInt();
    if(sql_query.value(1)!=NULL)
    {
        selectedTuple.case_dir          = sql_query.value(1).toString();
    }

    if(sql_query.value(2)!=NULL)
    {
        selectedTuple.case_create_date  = sql_query.value(2).toDateTime();
    }

    if(sql_query.value(3)!=NULL)
    {
        selectedTuple.case_operator    = sql_query.value(3).toString();
    }

    if(sql_query.value(4)!=NULL)
    {
        selectedTuple.patient_name      = sql_query.value(4).toString();
    }

    if(sql_query.value(5)!=NULL)
    {
        selectedTuple.patient_id        = sql_query.value(5).toString();
    }

    if(sql_query.value(6)!=NULL)
    {
        selectedTuple.patient_gender    = sql_query.value(6).toString();
    }

    if(sql_query.value(7)!=NULL)
    {
        selectedTuple.patient_age       = sql_query.value(7).toInt();
    }

    if(sql_query.value(8)!=NULL)
    {
        selectedTuple.surgery_operator  = sql_query.value(8).toString();
    }

    if (ok !=nullptr)
        *ok = true;

    return selectedTuple;
}

bool SqlDemo::SelectUser(QString pId, QString pPassword, user_info &pInfo)
{
    QSqlQuery q;
    q.prepare(SQL_SelectUserByIdAndPassword);
    q.addBindValue(pId);
    q.addBindValue(pPassword);
    if(!q.exec())
    {
        return false;
    }
    if(q.next())
    {
        int i=0;
        pInfo.id=q.value(i++).toString();
        pInfo.name=q.value(i++).toString();
        pInfo.password=q.value(i++).toString();
        pInfo.sex=q.value(i++).toInt();
        pInfo.authority=q.value(i++).toInt();
        pInfo.enable=q.value(i++).toInt();
        return true;
    }
    return false;
}

bool SqlDemo::IsUserExisted(QString pId,bool &pIsExisted)
{
    QSqlQuery q;
    q.prepare(SQL_IsUserExisted);
    q.addBindValue(pId);
    if(!q.exec())
    {
        return false;
    }
    pIsExisted=false;
    if(q.next())
        pIsExisted=true;
    return true;
}

bool SqlDemo::AddUser(const user_info &info)
{
    QSqlQuery q;
    q.prepare(SQL_InsertUser);
    q.addBindValue(info.id);
    q.addBindValue(info.name);
    q.addBindValue(encryptPassword(info.password));
    q.addBindValue(info.sex);
    q.addBindValue(info.authority);
    q.addBindValue(info.enable);
    if(!q.exec())
    {
        return false;
    }
    return true;
}

bool SqlDemo::GetUserByAuthority(std::vector<user_info> &pInfo, int pAuthority)
{
    QSqlQuery q;
    if(pAuthority==Authority_Admin)
    {
        q.prepare(SQL_SelectAllUser);
    }
    else
    {
        q.prepare(SQL_SelectUserByAuthority);
        q.addBindValue(pAuthority);
    }
    if(!q.exec())
    {
        return false;
    }
    pInfo.clear();
    while(q.next())
    {
        user_info info;
        int i=0;
        info.id=q.value(i++).toString();
        info.name=q.value(i++).toString();
        info.password=q.value(i++).toString();
        info.sex=q.value(i++).toInt();
        info.authority=q.value(i++).toInt();
        info.enable=q.value(i++).toInt();
        pInfo.push_back(info);
    }
    return true;
}

bool SqlDemo::GetUserByEnable(std::vector<user_info> &pInfo, Enable pEnable)
{
    QSqlQuery q;
    q.prepare(SQL_GetUserByEnable);
    q.addBindValue(pEnable);
    if(!q.exec())
        return false;
    while(q.next())
    {
        user_info info;
        int i=0;
        info.id=q.value(i++).toString();
        info.name=q.value(i++).toString();
        info.password=q.value(i++).toString();
        info.sex=q.value(i++).toInt();
        info.authority=q.value(i++).toInt();
        info.enable=q.value(i++).toInt();
        pInfo.push_back(info);
    }
    return true;
}

bool SqlDemo::RemoveUserById(QString pId)
{
    QSqlQuery q;
    q.prepare(SQL_RemoveUserById);
    q.addBindValue(pId);
    if(!q.exec())
    {
        return false;
    }
    return true;
}

bool SqlDemo::SetUserEnableById(QString pId, int pEnable)
{
    QSqlQuery q;
    q.prepare(SQL_SetUserEnableById);
    q.addBindValue(pEnable);
    q.addBindValue(pId);
    if(!q.exec())
        return false;
    return true;
}

bool SqlDemo::IsUserEnable(QString pId,int &pEnable)
{
    bool isExisted{false};
    if(!IsUserExisted(pId,isExisted) || !isExisted)
        return false;
    QSqlQuery q;
    q.prepare(SQL_IsUserEnable);
    q.addBindValue(pId);
    if(!q.exec())
        return false;
    pEnable=q.value(0).toInt();
    return true;
}

QMap<QString,QString> SqlDemo::SelectAllUserInfo()
{
    QString select_sql = QString("select id,password from usertable");
    QSqlQuery sql_query;
    sql_query.prepare(select_sql);
    QMap<QString,QString> userInfos;

    if (!sql_query.exec())
    {
        return userInfos;
    }

    while(sql_query.next())
    {
        QString userName=sql_query.value(0).toString();
        QString passWord=sql_query.value(1).toString();

        userInfos.insert(userName,passWord);
    }
    return userInfos;
}

std::vector<casetable_tuple> SqlDemo::SearchCaseByKeyword(QString keyword)
{
    QString strsql="select case_id,case_directory,case_create_date,patient_name,patient_id,patient_gender,patient_age,"
                     "surgery_operator,case_operator from casetable where case_id like '%keyword%' or patient_name like '%keyword%' order by case_id DESC;";
    strsql=strsql.replace("keyword",keyword);
    QSqlQuery sql_query;
    sql_query.prepare(strsql);
    std::vector<casetable_tuple> ret;
    if (!sql_query.exec())
    {
        return ret;
    }
    int i=0;
    casetable_tuple ca;
    while(sql_query.next())
    {
        i=0;
        ca.case_id=sql_query.value(i++).toInt();
        ca.case_dir=sql_query.value(i++).toInt();
        ca.case_create_date=sql_query.value(i++).toDateTime();
        ca.patient_name=sql_query.value(i++).toString();
        ca.patient_id=sql_query.value(i++).toString();
        ca.patient_gender=sql_query.value(i++).toString();
        ca.patient_age=sql_query.value(i++).toInt();
        ca.surgery_operator=sql_query.value(i++).toString();
        ca.case_operator=sql_query.value(i++).toString();
        ret.push_back(ca);
    }
    return ret;
}

std::vector<casetable_tuple> SqlDemo::SearchCaseByKeywordAndCaseOperator(QString keyword, QString pCaseOperator)
{
    QString strsql="select case_id,case_directory,case_create_date,patient_name,patient_id,patient_gender,patient_age,"
                     "surgery_operator from casetable where case_operator=? and (case_id like '%keyword%' or patient_name like '%keyword%' ) order by case_id DESC;";
    strsql=strsql.replace("keyword",keyword);
    QSqlQuery sql_query;
    sql_query.prepare(strsql);
    sql_query.addBindValue(pCaseOperator);
    std::vector<casetable_tuple> ret;
    if (!sql_query.exec())
    {
        //qDebug() << sql_query.lastError();
        return ret;
    }
    int i=0;
    casetable_tuple ca;
    while(sql_query.next())
    {
        i=0;
        ca.case_id=sql_query.value(i++).toInt();
        ca.case_dir=sql_query.value(i++).toInt();
        ca.case_create_date=sql_query.value(i++).toDateTime();
        ca.case_operator=pCaseOperator;
        ca.patient_name=sql_query.value(i++).toString();
        ca.patient_id=sql_query.value(i++).toString();
        ca.patient_gender=sql_query.value(i++).toString();
        ca.patient_age=sql_query.value(i++).toInt();
        ca.surgery_operator=sql_query.value(i++).toString();
        ret.push_back(ca);
    }
    return ret;
}

bool SqlDemo::RemoveCaseByCaseId(QString pCaseId)
{
    QString strsql=QString("delete from casetable where case_id = ? ");
    QSqlQuery sql_query;
    sql_query.prepare(strsql);
    sql_query.addBindValue(pCaseId);
    return sql_query.exec();
}


QString SqlDemo::encryptPassword(QString password)
{
    QString pwdmd5;
    QByteArray bb;
    bb = QCryptographicHash::hash(password.toLatin1(), QCryptographicHash::Md5);
    pwdmd5.append(bb.toHex());
    return pwdmd5;
}

bool SqlDemo::UpdataLoginTableTuple(QString id, QString password)
{
    QString pwdmd5=encryptPassword(password);

    QString updata_sql = QString("update usertable set password = ? where id = ?");
    QSqlQuery sql_query;
    sql_query.prepare(updata_sql);
    sql_query.addBindValue(password);
    sql_query.addBindValue(id);
    if (!sql_query.exec())
    {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;

}

bool SqlDemo::UpdataLoginTableTuple(QString id, QString password, bool encrypt)
{
    QString pwd=password;
    if(encrypt)
        pwd=encryptPassword(password);
    QString updata_sql;
    if(pwd.isEmpty())
    {
        updata_sql = QString("update usertable set password = NULL where id = '%1'").arg(id);
    }
    else
    {
        updata_sql = QString("update usertable set password = '%1' where id = '%2'").arg(pwd).arg(id);
    }
    QSqlQuery sql_query;
    sql_query.prepare(updata_sql);
    if (!sql_query.exec())
    {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
