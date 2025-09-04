/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#ifdef SUANPAN_VTK

// ReSharper disable StringLiteralTypo
#include "vtkParser.h"

#include <Domain/DomainBase.h>
#include <Domain/Node.h>
#include <Element/Element.h>
#include <Toolbox/utility.h>
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2)  // NOLINT(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)
VTK_MODULE_INIT(vtkInteractionStyle)  // NOLINT(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)
VTK_MODULE_INIT(vtkRenderingFreeType) // NOLINT(cppcoreguidelines-special-member-functions, hicpp-special-member-functions)

#include <vtkActor.h>
#include <vtkColorTransferFunction.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkLookupTable.h>
#include <vtkNamedColors.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkScalarBarActor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>

vtkInfo vtk_process(std::istringstream& command) {
    vtkInfo config;

    std::string keyword;

    while(!command.eof() && get_input(command, keyword))
        if(is_equal(keyword, "scale") && !get_input(command, config.scale)) config.scale = 1.;
        else if(is_equal(keyword, "deformed")) config.on_deformed = true;
        else if(is_equal(keyword, "undeformed")) config.on_deformed = false;
        else if(is_equal(keyword, "type") && get_input(command, keyword)) config.type = to_token(keyword);
        else if(is_equal(keyword, "fontsize") && !get_input(command, config.font_size)) config.font_size = 8;
        else if(is_equal(keyword, "save") && get_input(command, config.file_name)) config.save_file = true;
        else if(is_equal(keyword, "nobar")) config.colorbar = false;
        else if(is_equal(keyword, "material") && !get_input(command, config.material_type)) config.material_type = -1;
        else if(is_equal(keyword, "size")) {
            if(!get_input(command, config.canvas_size[0])) config.canvas_size[0] = 500;
            if(!get_input(command, config.canvas_size[1])) config.canvas_size[1] = 500;
        }

    return config;
}

void vtk_setup(const vtkSmartPointer<vtkUnstructuredGrid>& grid, const vtkInfo& config) {
    const auto color = vtkSmartPointer<vtkNamedColors>::New();
    const auto table = vtkSmartPointer<vtkLookupTable>::New();
    const auto func = vtkSmartPointer<vtkColorTransferFunction>::New();
    const auto mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    const auto bar = vtkSmartPointer<vtkScalarBarActor>::New();
    const auto actor = vtkSmartPointer<vtkActor>::New();
    const auto interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    const auto renderer = vtkSmartPointer<vtkRenderer>::New();
    const auto window = vtkSmartPointer<vtkRenderWindow>::New();

    func->SetColorSpaceToDiverging();
    func->AddRGBPoint(0., .230, .299, .754);
    func->AddRGBPoint(1., .706, .016, .150);
    table->SetNumberOfTableValues(256);
    table->Build();
    double rgb[4] = {0., 0., 0., 1.};
    for(auto I = 0; I < 256; ++I) {
        func->GetColor(static_cast<double>(I) / 256., rgb);
        table->SetTableValue(I, rgb);
    }

    mapper->SetInputDataObject(grid);
    mapper->SetLookupTable(table);
    mapper->SetScalarRange(grid->GetPointData()->GetScalars()->GetRange());

    bar->SetLookupTable(mapper->GetLookupTable());

    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(color->GetColor3d("DodgerBlue").GetData());
    actor->GetProperty()->EdgeVisibilityOn();
    actor->GetProperty()->SetOpacity(1.);

    renderer->AddActor(actor);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    if(config.colorbar) renderer->AddActor2D(bar);
#pragma GCC diagnostic pop
    renderer->SetBackground(color->GetColor3d("Grey").GetData());
    renderer->ResetCameraClippingRange();

    window->AddRenderer(renderer);

    interactor->SetRenderWindow(window);

    window->Render();
    window->SetSize(config.canvas_size[0], config.canvas_size[1]);
    window->SetWindowName(config.title_name.c_str());

    interactor->Start();
}

void vtk_save(vtkSmartPointer<vtkUnstructuredGrid>&& grid, const vtkInfo config) {
    // NOLINT(performance-unnecessary-value-param)
    const auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetInputData(grid);
    writer->SetFileName(config.file_name.c_str());
    writer->SetFileTypeToBinary();
    writer->Write();
    suanpan_debug("Plot is written to file \"{}\".\n", config.file_name);
}

int vtk_parser(const shared_ptr<DomainBase>& domain, std::istringstream& command) {
    auto plot_info = vtk_process(command);
    if(!plot_info.on_deformed) plot_info.scale = 0.;

    const auto L = to_token(to_category(plot_info.type));

    const auto func = OutputType::U == L || OutputType::V == L || OutputType::A == L || OutputType::RF == L || OutputType::DF == L || OutputType::IF == L ? vtk_plot_node_quantity : vtk_plot_element_quantity;

#ifdef SUANPAN_WIN
    domain->insert(std::async(std::launch::async, func, std::cref(domain), plot_info));
#else
    func(domain, plot_info);
#endif

    return SUANPAN_SUCCESS;
}

void vtk_plot_node_quantity(const shared_ptr<DomainBase>& domain, vtkInfo config) {
    auto& t_node_pool = domain->get_node_pool();
    auto& t_element_pool = domain->get_element_pool();

    config.title_name = "Plotting Nodal Quantity " + std::string(to_name(config.type));

    auto max_node = static_cast<unsigned>(t_node_pool.size());
    for(const auto& I : t_node_pool) max_node = std::max(max_node, I->get_tag());

    auto data = vtkSmartPointer<vtkDoubleArray>::New();
    auto node = vtkSmartPointer<vtkPoints>::New();

    data->SetNumberOfComponents(6);
    data->SetNumberOfTuples(++max_node);
    node->SetNumberOfPoints(max_node);

    for(unsigned I = 0; I < max_node; ++I) {
        node->SetPoint(I, 0., 0., 0.);
        data->SetTuple6(I, 0., 0., 0., 0., 0., 0.);
    }

    data->SetName(to_category(config.type).c_str());
    data->SetComponentName(0, "1");
    data->SetComponentName(1, "2");
    data->SetComponentName(2, "3");
    data->SetComponentName(3, "4");
    data->SetComponentName(4, "5");
    data->SetComponentName(5, "6");

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->Allocate(static_cast<vtkIdType>(t_element_pool.size()));
    std::ranges::for_each(t_element_pool, [&](const shared_ptr<Element>& t_element) {
        t_element->SetDeformation(node, config.scale);
        t_element->GetData(data, to_token(to_category(config.type)));
        if(const auto t_cell = t_element->GetCell(); nullptr != t_cell) grid->InsertNextCell(t_cell->GetCellType(), t_cell->GetPointIds());
    });

    grid->SetPoints(node);

    if(config.store_ptr) {
        grid->GetPointData()->SetScalars(data);
        grid->GetPointData()->SetActiveScalars(to_category(config.type).c_str());
        config.grid_ptr = grid;
    }
    else if(config.save_file) {
        grid->GetPointData()->SetScalars(data);
        grid->GetPointData()->SetActiveScalars(to_category(config.type).c_str());
        domain->insert(std::async(std::launch::async, vtk_save, std::move(grid), std::move(config)));
    }
    else {
        const auto sub_data = vtkSmartPointer<vtkDoubleArray>::New();

        sub_data->SetNumberOfTuples(data->GetNumberOfTuples());
        sub_data->CopyComponent(0, data, to_index(config.type));

        grid->GetPointData()->SetScalars(sub_data);

        vtk_setup(grid, config);
    }
}

void vtk_plot_element_quantity(const shared_ptr<DomainBase>& domain, vtkInfo config) {
    auto& t_node_pool = domain->get_node_pool();
    auto& t_element_pool = domain->get_element_pool();

    config.title_name = "Plotting Element Quantity " + std::string(to_name(config.type));

    auto max_node = static_cast<unsigned>(t_node_pool.size());
    for(const auto& I : t_node_pool) max_node = std::max(max_node, I->get_tag());

    const auto data = vtkSmartPointer<vtkDoubleArray>::New();
    auto node = vtkSmartPointer<vtkPoints>::New();

    data->SetNumberOfComponents(6);
    data->SetNumberOfTuples(++max_node);
    node->SetNumberOfPoints(max_node);

    for(unsigned I = 0; I < max_node; ++I) {
        node->SetPoint(I, 0., 0., 0.);
        data->SetTuple6(I, 0., 0., 0., 0., 0., 0.);
    }

    data->SetName(to_category(config.type).c_str());
    data->SetComponentName(0, "1");
    data->SetComponentName(1, "2");
    data->SetComponentName(2, "3");
    data->SetComponentName(3, "4");
    data->SetComponentName(4, "5");
    data->SetComponentName(5, "6");

    mat tensor(6, max_node, fill::zeros);
    vec counter(max_node, fill::zeros);

    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    grid->Allocate(static_cast<vtkIdType>(t_element_pool.size()));
    std::ranges::for_each(t_element_pool, [&](const shared_ptr<Element>& t_element) {
        if(-1 != config.material_type)
            if(const auto& t_tag = t_element->get_material_tag(); t_tag.empty() || static_cast<uword>(config.material_type) != t_tag(0)) return;
        t_element->SetDeformation(node, config.scale);
        auto& t_encoding = t_element->get_node_encoding();
        counter(t_encoding) += 1.;
        if(const auto t_data = t_element->GetData(to_token(to_category(config.type))); !t_data.empty()) tensor.cols(t_encoding) += t_data;
        if(const auto t_cell = t_element->GetCell(); nullptr != t_cell) grid->InsertNextCell(t_cell->GetCellType(), t_cell->GetPointIds());
    });

    for(unsigned I = 0; I < max_node; ++I)
        if(0. != counter(I)) {
            tensor.col(I) /= counter(I);
            data->SetTuple(I, tensor.colptr(I));
        }

    grid->SetPoints(node);

    if(config.store_ptr) {
        grid->GetPointData()->SetScalars(data);
        grid->GetPointData()->SetActiveScalars(to_category(config.type).c_str());
        config.grid_ptr = grid;
    }
    else if(config.save_file) {
        grid->GetPointData()->SetScalars(data);
        grid->GetPointData()->SetActiveScalars(to_category(config.type).c_str());
        domain->insert(std::async(std::launch::async, vtk_save, std::move(grid), std::move(config)));
    }
    else {
        const auto sub_data = vtkSmartPointer<vtkDoubleArray>::New();

        sub_data->SetNumberOfTuples(data->GetNumberOfTuples());
        sub_data->CopyComponent(0, data, to_index(config.type));

        grid->GetPointData()->SetScalars(sub_data);

        vtk_setup(grid, config);
    }
}

#endif
