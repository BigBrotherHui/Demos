#ifndef STRSQL_H
#define STRSQL_H
#include <QString>
//usertable
static const QString SQL_CreateLoginTable="create table if not exists usertable (id varchar2(200) primary key,name varchar2(200),password varchar2(200),"
                                           "sex int,authority INT DEFAULT (0),enable int default (0))";
static const QString SQL_CreateCaseTable = "create table if not exists casetable (case_id integer primary key autoincrement, case_directory varchar2(200), \
                    case_create_date varchar2(200),case_operator varchar2(200),patient_name varchar2(30), patient_id varchar2(30), \
                    patient_gender varchar2(1), patient_age int, surgery_operator VARCHAR2(30), \
                    foreign key (surgery_operator) references usertable(id))";
static const QString SQL_InsertUser = "insert into usertable values (?,?,?,?,?,?)";
static const QString SQL_SelectUserByIdAndPassword ="select * from usertable where id=? and password=?";
static const QString SQL_IsUserExisted="select id from usertable where id = ?";
static const QString SQL_SelectUserByAuthority="select * from usertable where authority=?";
static const QString SQL_SelectAllUser="select * from usertable";
static const QString SQL_RemoveUserById="delete from usertable where id=?";
static const QString SQL_SetUserEnableById="update usertable set enable=? where id=?";
static const QString SQL_IsUserEnable="select enable from usertable where id = ?";
static const QString SQL_UpdateUserPasswordById="update usertable set password=? where id=?";
static const QString SQL_GetUserByEnable="select * from usertable where enable=?";














#endif // STRSQL_H
