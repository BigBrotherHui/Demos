#pragma execution_character_set("utf-8")
#include "cryptoutil.h"
#include "stdio.h"
#include "QDebug"

#include <QCryptographicHash>
#include <QTextCodec>

namespace CryptoUtil{
const int ALPHABET_SIZE = 64;
const char ALPHABET[] = "ABCDEFGH" "IJKLMNOP" "QRSTUVWX" "YZabcdef"
                               "ghijklmn" "opqrstuv" "wxyz0123" "456789+/";
const char PADCHAR = '=';
const char PADCHAR_REVERSE = '~';
//base64加密
QString base64Encrypt(QString code){
    int padlen = 0;
    std::string sourceData = code.toStdString();
    std::string tmp;
    tmp.resize((sourceData.size() + 2) / 3 * 4);

    int i = 0;
    char *out = &tmp[0];
    while (i < sourceData.size()) {
        // encode 3 bytes at a time
        int chunk = 0;
        chunk |= int(uchar(sourceData.data()[i++])) << 16;
        if (i == sourceData.size()) {
            padlen = 2;
        } else {
            chunk |= int(uchar(sourceData.data()[i++])) << 8;
            if (i == sourceData.size())
                padlen = 1;
            else
                chunk |= int(uchar(sourceData.data()[i++]));
        }

        int j = (chunk & 0x00fc0000) >> 18;
        int k = (chunk & 0x0003f000) >> 12;
        int l = (chunk & 0x00000fc0) >> 6;
        int m = (chunk & 0x0000003f);
        *out++ = ALPHABET[j];
        *out++ = ALPHABET[k];

        if (padlen > 1) {
            *out++ = PADCHAR;
        } else {
            *out++ = ALPHABET[l];
        }
        if (padlen > 0) {
            *out++ = PADCHAR;
        } else {
            *out++ = ALPHABET[m];
        }
    }
    assert(/*(options & OmitTrailingEquals)*/0 || (out == tmp.size() + tmp.data()));
    if (/*options & OmitTrailingEquals*/0)
        tmp.resize(out - tmp.data());
    return QString::fromStdString(tmp);
}

//base64解密
QString base64Decrypt(QString code){

    unsigned int buf = 0;
    int nbits = 0;
    std::string tmp;
    std::string sourceData = code.toStdString();
    tmp.resize((sourceData.size() * 3) / 4);

    int offset = 0;
    for (int i = 0; i < sourceData.size(); ++i) {
        int ch = sourceData.at(i);
        int d;

        if (ch >= 'A' && ch <= 'Z')
            d = ch - 'A';
        else if (ch >= 'a' && ch <= 'z')
            d = ch - 'a' + 26;
        else if (ch >= '0' && ch <= '9')
            d = ch - '0' + 52;
        else if (ch == '+')
            d = 62;
        else if (ch == '/')
            d = 63;
        else
            d = -1;

        if (d != -1) {
            buf = (buf << 6) | d;
            nbits += 6;
            if (nbits >= 8) {
                nbits -= 8;
                tmp[offset++] = buf >> nbits;
                buf &= (1 << nbits) - 1;
            }
        }
    }

    tmp.resize(offset);
    return QString::fromStdString(tmp);
}

// 混淆
QString confusion(QString code){
    //将字符串从中间分割为2部分，将后半部分继续进行混淆后，按顺序逐个插入前半部分中
    int len = code.length();
    if(len<=1){
        return code;
    }
    int leftLen = (len + 1) / 2;
    QString right = code.mid(leftLen);
    QString newRight = confusion(right);
    QString newCode;
    for(int i=0;i<leftLen;i++){
        newCode += code.mid(i,1);
        if(i < newRight.length()){
            newCode += newRight.mid(i,1);
        }
    }
    return newCode;
}

// 反混淆
QString antiConfusion(QString code){
    //字符串中奇数位为最终结果的前半部分，偶数位为混淆后的最终结果后半部分
    int len = code.length();
    if(len<=1){
        return code;
    }
    QString left, right;
    for(int i=0;i<len;i++){
        left += code.mid(i,1);
        if(i < len - 1){
            right += code.mid(++i,1);
        }
    }
    QString newCode = left + antiConfusion(right);
    return newCode;
}

char reverseChar(char ch){
    for(int j=0;j<ALPHABET_SIZE;j++){
        if(ch==ALPHABET[j]){
            return ALPHABET[ALPHABET_SIZE - j - 1];
        }
    }
    return ch;
}
// 反转
QString reverse(QString code){
    if(code.length()<=0){
        return code;
    }
    QByteArray array = code.toUtf8();
    QByteArray newArray;
    for(int i=0;i<array.size();i++){
        char ch = array.at(i);
        if(ch==PADCHAR){
            newArray[i] = PADCHAR_REVERSE;
        }
        else if(ch==PADCHAR_REVERSE){
            newArray[i] = PADCHAR;
        }
        else{
            newArray[i] = reverseChar(ch);
        }
    }
    QString newCode = newArray;
    return newCode;
}

// 加密
QString encrypt(QString code){
    //对明码进行混淆后进行BASE64加密，然后反转
    return reverse(base64Encrypt(confusion(code)));
}

// 解密
QString decrypt(QString code){
    //对密码进行反转后，使用base64解密后解除混淆
    return antiConfusion(base64Decrypt(reverse(code)));
}

// md5加密
QString md5(QString code){
    return QCryptographicHash::hash(code.toUtf8(),QCryptographicHash::Md5).toHex();
}
}
