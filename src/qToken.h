// Author: Jerry Chou
//      Lawrence Berkeley National Laboratory
#ifndef FQ_QTOKEN_H
#define FQ_QTOKEN_H
#include <string>

class QToken{
public:
    enum TYPE {
	VAR, NUMBER, ARTHIMETICOP, RANGEOP, LOGICOP
    };

    QToken(TYPE type, const char* text) {
	this->type = type;
	this->text = text;
	this->left = 0;
	this->right = 0;
    }
    ~QToken() {delete right; delete left;}
    TYPE getType() {return type;}
    void setLeft(QToken* qToken) {left = qToken;}
    void setRight(QToken* qToken) {right = qToken;}
    QToken* getLeft() {return left;}
    QToken* getRight() {return right;}
    std::string getText(){return text;}
private:
    TYPE type;
    std::string text;
    QToken* left;
    QToken* right;
};
#endif
