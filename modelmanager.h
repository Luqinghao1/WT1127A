#ifndef MODELMANAGER_H
#define MODELMANAGER_H

#include <QObject>
#include <QWidget>
#include <QStackedWidget>
#include <QComboBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QMap>
#include <QVector>
#include <tuple>

class ModelWidget1;
class ModelWidget2;
class ModelWidget3;

// 定义计算结果类型：Time, Pressure, Derivative
typedef std::tuple<QVector<double>, QVector<double>, QVector<double>> ModelCurveData;

class ModelManager : public QObject
{
    Q_OBJECT

public:
    enum ModelType {
        InfiniteConductive = 0,    // 无限导流双重孔隙介质页岩油藏
        FiniteConductive = 1,      // 有限导流 (预留)
        SegmentedMultiCluster = 2  // 分段多簇 (预留)
    };
    Q_ENUM(ModelType)

    explicit ModelManager(QWidget* parent = nullptr);
    ~ModelManager();

    void initializeModels(QWidget* parentWidget);

    QWidget* getMainWidget() const { return m_mainWidget; }
    QStackedWidget* getModelStack() const { return m_modelStack; }
    void switchToModel(ModelType modelType);
    ModelType getCurrentModelType() const { return m_currentModelType; }

    static QString getModelTypeName(ModelType modelType);
    static QStringList getAvailableModelTypes();

    ModelWidget1* getInfiniteConductive() const { return m_modelWidget1; }
    ModelWidget2* getFiniteConductive() const { return m_modelWidget2; }
    ModelWidget3* getSegmentedMultiCluster() const { return m_modelWidget3; }

    // --- 计算接口 ---
    QMap<QString, double> getDefaultParameters(ModelType type);
    ModelCurveData calculateTheoreticalCurve(ModelType type, const QMap<QString, double>& params);

signals:
    void modelSwitched(ModelType newModelType, ModelType oldModelType);
    void calculationCompleted(const QString &analysisType, const QMap<QString, double> &results);

private slots:
    void onModelTypeSelectionChanged(int index);
    void onModel1CalculationCompleted(const QString &analysisType, const QMap<QString, double> &results);
    void onModel2CalculationCompleted(const QString &analysisType, const QMap<QString, double> &results);
    void onModel3CalculationCompleted(const QString &analysisType, const QMap<QString, double> &results);

private:
    void setupModelSelection();
    void connectModelSignals();
    void createMainWidget();

    // 具体模型入口
    ModelCurveData calculateModel1(const QMap<QString, double>& params);
    ModelCurveData calculateModel2(const QMap<QString, double>& params);
    ModelCurveData calculateModel3(const QMap<QString, double>& params);

    // --- Stehfest 算法数学核心 (移植自 ModelWidget1，确保曲线一致) ---
    QVector<double> solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b);
    double stefestCoefficient(int i, int N);
    double factorial(int n);

    // 拉普拉斯空间解
    double flaplace(double z, const QMap<QString, double>& params);
    double e_function(double z, int i, int j, int k, int v, int mf, int nf,
                      double omega, double lambda, double Xf, double yy, double y);
    double f_function(int j, int nf, double Xf, double y);

    // 贝塞尔函数与积分
    double integralBesselK0(double XDkv, double YDkv, double yDij, double fz, double xDij, double xDij1);
    double besselK0(double x);
    double gaussQuadrature(double XDkv, double YDkv, double yDij, double fz, double a, double b);

    // 压力导数计算
    void computePressureDerivative(const QVector<double>& tD, const QVector<double>& pd,
                                   double cD, QVector<double>& dpd, QVector<double>& td_dpd);

private:
    QWidget* m_mainWidget;
    QComboBox* m_modelTypeCombo;
    QStackedWidget* m_modelStack;

    ModelWidget1* m_modelWidget1;
    ModelWidget2* m_modelWidget2;
    ModelWidget3* m_modelWidget3;

    ModelType m_currentModelType;
};

#endif // MODELMANAGER_H
