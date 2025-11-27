#ifndef FITTINGWIDGET_H
#define FITTINGWIDGET_H

#include <QWidget>
#include <QMap>
#include <QVector>
#include <QFuture>
#include <QFutureWatcher>
#include <QDialog>
#include <QComboBox>
#include <QTableWidget>
#include <QLabel>
#include <QPushButton>
#include "modelmanager.h"
#include "qcustomplot.h"

namespace Ui { class FittingWidget; }

// ===========================================================================
// 数据加载对话框类
// ===========================================================================
class FittingDataLoadDialog : public QDialog
{
    Q_OBJECT
public:
    explicit FittingDataLoadDialog(const QList<QStringList>& previewData, QWidget *parent = nullptr);
    int getTimeColumnIndex() const;
    int getPressureColumnIndex() const;
    int getDerivativeColumnIndex() const;
    int getSkipRows() const;
private slots:
    void validateSelection();
private:
    QTableWidget* m_previewTable;
    QComboBox* m_comboTime, *m_comboPressure, *m_comboDeriv, *m_comboSkipRows;
};

// ===========================================================================
// 拟合参数结构体
// ===========================================================================
struct FitParameter {
    QString name;       // 内部参数名 (如 "omega")
    QString displayName;// 显示名称 (如 "储容比 (ω)")
    double value;       // 当前值
    double min;         // 最小值下限
    double max;         // 最大值上限
    bool isFit;         // 是否参与拟合
};

// ===========================================================================
// 拟合主窗口类
// ===========================================================================
class FittingWidget : public QWidget
{
    Q_OBJECT
public:
    explicit FittingWidget(QWidget *parent = nullptr);
    ~FittingWidget();

    // 设置模型管理器
    void setModelManager(ModelManager* manager);
    // 设置观测数据（实测数据）
    void setObservedData(const QVector<double>& time, const QVector<double>& pressure, const QVector<double>& derivative);

signals:
    // 拟合完成信号
    void fittingCompleted(ModelManager::ModelType modelType, const QMap<QString, double> &parameters);
    // 迭代更新信号（用于刷新界面上的误差和曲线）
    void sigIterationUpdated(double currentError, const QMap<QString, double>& currentParams);
    // 进度信号
    void sigProgress(int value);

private slots:
    // 模型选择改变时触发
    void on_comboModelSelect_currentIndexChanged(int index);
    // 重置参数按钮点击
    void on_btnResetParams_clicked();
    // 开始拟合按钮点击
    void on_btnRunFit_clicked();
    // 停止按钮点击
    void on_btnStop_clicked();
    // 应用结果按钮点击
    void on_btnUpdateModel_clicked();
    // 加载数据按钮点击
    void on_btnLoadData_clicked();

    // 跨线程更新界面的槽函数
    void onIterationUpdate(double currentError, const QMap<QString, double>& currentParams);
    // 拟合任务结束时的槽函数
    void onFitFinished();

private:
    Ui::FittingWidget *ui;
    QCustomPlot *m_plot;
    ModelManager* m_modelManager;

    // 观测数据（实测）
    QVector<double> m_obsTime, m_obsPressure, m_obsDerivative;
    // 参数列表
    QList<FitParameter> m_parameters;

    bool m_isFitting;     // 是否正在拟合中
    bool m_stopRequested; // 是否请求停止
    QFutureWatcher<void> m_watcher; // 异步任务监视器

    // --- 界面初始化与辅助函数 ---
    void setupPlot();       // 初始化绘图控件
    void initModelCombo();  // 初始化模型下拉框
    void loadParamsToTable(); // 将参数加载到表格显示
    void updateParamsFromTable(); // 从表格读取参数回内存

    // 绘制曲线
    void plotCurves(const QVector<double>& t, const QVector<double>& p, const QVector<double>& dp, bool isModel);

    // 获取参数的中文显示名称
    QString getParamDisplayName(const QString& key);

    // --- 核心优化逻辑 ---
    // [修复] 增加参数传递，避免子线程直接访问 UI
    void runOptimizationTask(int algIndex, ModelManager::ModelType modelType);

    // 具体算法实现
    void runLevenbergMarquardtOptimization(ModelManager::ModelType modelType); // LM 算法
    void runNelderMeadOptimization(ModelManager::ModelType modelType);         // 单纯形法

    // 数学辅助函数
    // [修复] calculateResiduals 需要传入 modelType，不能从 UI 获取
    QVector<double> calculateResiduals(const QMap<QString, double>& params, ModelManager::ModelType modelType);
    double calculateSumSquaredError(const QVector<double>& residuals);

    // 计算雅可比矩阵
    QVector<QVector<double>> computeJacobian(const QMap<QString, double>& params, const QVector<double>& currentResiduals, const QVector<int>& fitIndices, ModelManager::ModelType modelType);
    // 求解线性方程组
    QVector<double> solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b);

    QStringList parseLine(const QString& line); // 解析文件行
    // [修复] calculateError 也需要传入 modelType
    double calculateError(const QMap<QString, double>& trialParams, ModelManager::ModelType modelType);
};
#endif // FITTINGWIDGET_H
