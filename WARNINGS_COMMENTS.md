# Комментарии к Warnings в тестах

## Анализ warnings из `test_finder.py`

### 1. DeprecationWarning: pkg_resources is deprecated (seqfold)
```
../../../../../biograph/prod/.venv/lib/python3.10/site-packages/seqfold/__init__.py:2
DeprecationWarning: pkg_resources is deprecated as an API
```

**Причина:** Внешняя зависимость `seqfold` использует устаревший API `pkg_resources` для получения версии пакета.

**Статус:** ⚠️ Не критично, но стоит отслеживать
- Это предупреждение от внешней библиотеки, не от нашего кода
- `seqfold` нужно обновить до версии, использующей современный API
- Не влияет на функциональность, но может исчезнуть в будущих версиях Python

**Рекомендация:** 
- Отслеживать обновления `seqfold`
- При обновлении зависимостей проверить, исправлена ли проблема

---

### 2. DeprecationWarning: Deprecated call to `pkg_resources.declare_namespace` (pkg_resources)
```
../../../../../biograph/prod/.venv/lib/python3.10/site-packages/pkg_resources/__init__.py:3142
DeprecationWarning: Deprecated call to `pkg_resources.declare_namespace('sphinxcontrib')`
```

**Причина:** Другая внешняя зависимость (возможно `sphinxcontrib` или связанная) использует устаревший метод объявления namespace-пакетов.

**Статус:** ⚠️ Не критично
- Это предупреждение от внешней зависимости, не от нашего кода
- Связано с namespace packages в Python
- Не влияет на функциональность

**Рекомендация:**
- Игнорировать до обновления зависимостей
- При обновлении зависимостей проверить, исправлена ли проблема

---

### 3. UserWarning: Hairpin energy threshold (наше предупреждение)
```
/biograph/prod/.venv/lib/python3.10/site-packages/ollikit/oligama/exceptions.py:30
UserWarning: Hairpin energy threshold (0.00 kJ/mol) is greater than hairpin energy for s1 (-3.70 kJ/mol). Algorithm may not converge.
```

**Причина:** Это **наше** предупреждение, которое мы намеренно генерируем в коде.

**Статус:** ✅ Ожидаемое поведение
- Это информационное предупреждение для пользователя
- Указывает, что порог энергии шпильки может быть слишком высоким
- Алгоритм может не сойтись или работать медленнее

**Рекомендация:**
- Это нормальное поведение для теста
- В тесте используется `Hairpin_energy_thr: 0`, что может быть слишком строгим условием
- Можно оставить как есть или настроить более реалистичные значения в тесте

**Код предупреждения:** `ollikit/oligama/scripts/olig_finder.py:15-19`

---

### 4. UserWarning: A worker stopped (joblib)
```
/biograph/prod/.venv/lib/python3.10/site-packages/joblib/externals/loky/process_executor.py:752
UserWarning: A worker stopped while some jobs were given to the executor. This can be caused by a too short worker timeout or by a memory leak.
```

**Причина:** Проблема с параллельной обработкой через `joblib` (используется в предикторах Nupack).

**Статус:** ⚠️ Требует внимания
- Может указывать на проблему с таймаутами воркеров
- Или на утечку памяти в параллельных процессах
- Или на слишком короткий таймаут для воркеров

**Возможные причины:**
1. **Слишком короткий таймаут воркера** - задачи выполняются дольше, чем ожидается
2. **Утечка памяти** - воркеры потребляют слишком много памяти и завершаются
3. **Проблемы с ресурсами** - нехватка CPU/памяти для параллельной обработки

**Рекомендация:**
- Проверить настройки `joblib` в коде предикторов
- Увеличить таймаут воркеров, если проблема в таймауте
- Проверить использование памяти при параллельной обработке
- Возможно, стоит уменьшить количество параллельных воркеров

**Где искать:**
- `ollikit/oligama/predictors/` - предикторы, использующие Nupack
- Настройки `Parallel` в `joblib`

---

## Итоговая оценка

| Warning | Критичность | Действие |
|---------|-------------|----------|
| seqfold pkg_resources | Низкая | Отслеживать обновления |
| pkg_resources namespace | Низкая | Игнорировать |
| Hairpin energy threshold | Нет (ожидаемое) | Оставить как есть |
| joblib worker stopped | Средняя | Проверить настройки joblib |

## Рекомендации по подавлению warnings в тестах

Если нужно подавить некоторые warnings в тестах, можно использовать:

```python
import warnings
import pytest

# Подавить deprecation warnings от внешних библиотек
@pytest.fixture(autouse=True)
def suppress_external_warnings():
    warnings.filterwarnings("ignore", category=DeprecationWarning, module="seqfold")
    warnings.filterwarnings("ignore", category=DeprecationWarning, module="pkg_resources")
    # Но НЕ подавлять наши UserWarning - они информативны
```

Или в `pytest.ini` или `pyproject.toml`:
```toml
[tool.pytest.ini_options]
filterwarnings = [
    "ignore::DeprecationWarning:seqfold.*",
    "ignore::DeprecationWarning:pkg_resources.*",
]
```

