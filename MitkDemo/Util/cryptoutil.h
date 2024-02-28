#ifndef CRYPTOUTIL_H
#define CRYPTOUTIL_H

#include <QString>

namespace CryptoUtil
{
    /**
     * @brief 加密处理
     * @param code 加密前的明文
     * @return 加密后的密文
     */
    QString encrypt(QString code);
    /**
     * @brief 解密处理
     * @param code 解密前的密文
     * @return 解密后的明文
     */
    QString decrypt(QString code);
    /**
     * @brief MD5加密处理
     * @param code 加密前的明文
     * @return 加密后的密文
     */
    QString md5(QString code);
    /**
     * @brief BASE64加密处理
     * @param code 加密前的明文
     * @return 加密后的密文
     */
    QString base64Encrypt(QString code);
    /**
     * @brief BASE64解密处理
     * @param code 解密前的密文
     * @return 解密后的明文
     */
    QString base64Decrypt(QString code);
};

#endif // CRYPTOUTIL_H
