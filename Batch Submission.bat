@echo off

REM Submit the First job
echo Wait for the First job to complete, do not close this window...
start /wait cmd /C abaqus job=One_Element user=One_Element int ask=off

REM Submit the Second job
echo Wait for the Second job to complete, do not close this window...
start /wait cmd /C abaqus job=One_Element user=One_Element int ask=off

REM Submit the Third job
echo Wait for the Third job to complete, do not close this window...
start /wait cmd /C abaqus job=One_Element user=One_Element int ask=off


