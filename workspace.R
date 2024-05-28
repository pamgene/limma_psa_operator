

rm(list = ls())
library(tercen)
library(dplyr)

options("tercen.workflowId" = "d81a58f796a0528726a004cfc30294b8")
options("tercen.stepId"     = "3aa0a5dd-531c-46d4-9f9b-f6fdc7c94a63")

getOption("tercen.workflowId")
getOption("tercen.stepId")

source("main.R")
